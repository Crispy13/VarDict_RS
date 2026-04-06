"""Reporting helpers for parity sweeps."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import json
import xml.etree.ElementTree as ET


@dataclass(frozen=True, slots=True)
class ShardFailure:
    shard_label: str
    region: str
    reason: str
    first_diff_line: int | None = None
    summary: str = ""


@dataclass(slots=True)
class ChromResult:
    status: str
    total_shards: int
    pass_count: int
    fail_count: int
    empty_count: int
    failures: list[ShardFailure] = field(default_factory=list)


@dataclass(slots=True)
class ConfigResult:
    config_id: str
    label: str
    tier: int | str
    chromosomes: dict[str, ChromResult]
    elapsed_seconds: float


def print_summary_table(results: list[ConfigResult]) -> None:
    headers = ("Config", "Label", "Tier", "Chroms", "Pass", "Fail", "Empty", "Elapsed", "Status")
    rows: list[tuple[str, ...]] = []
    for result in results:
        chrom_count = len(result.chromosomes)
        pass_count = sum(chrom.pass_count for chrom in result.chromosomes.values())
        fail_count = sum(chrom.fail_count for chrom in result.chromosomes.values())
        empty_count = sum(chrom.empty_count for chrom in result.chromosomes.values())
        rows.append(
            (
                result.config_id,
                result.label,
                str(result.tier),
                str(chrom_count),
                str(pass_count),
                str(fail_count),
                str(empty_count),
                f"{result.elapsed_seconds:.1f}s",
                _config_status(result),
            )
        )

    widths = [len(header) for header in headers]
    for row in rows:
        widths = [max(width, len(value)) for width, value in zip(widths, row, strict=True)]

    print(_format_row(headers, widths))
    print(_format_row(tuple("-" * width for width in widths), widths))
    for row in rows:
        print(_format_row(row, widths))


def print_status_matrix(status: dict[str, dict[str, str]]) -> None:
    chromosomes = sorted({chrom for row in status.values() for chrom in row})
    headers = ("Config", *chromosomes)
    rows: list[tuple[str, ...]] = []
    for config_id in sorted(status):
        row = [config_id]
        for chrom in chromosomes:
            row.append(status[config_id].get(chrom, "PEND"))
        rows.append(tuple(row))

    widths = [len(header) for header in headers]
    for row in rows:
        widths = [max(width, len(value)) for width, value in zip(widths, row, strict=True)]

    print(_format_row(headers, widths))
    print(_format_row(tuple("-" * width for width in widths), widths))
    for row in rows:
        print(_format_row(row, widths))


def write_json_report(results: list[ConfigResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "results": [_config_result_to_dict(result) for result in results],
        "summary": {
            "configs": len(results),
            "pass": sum(1 for result in results if _config_status(result) == "PASS"),
            "fail": sum(1 for result in results if _config_status(result) == "FAIL"),
            "pending": sum(1 for result in results if _config_status(result) in {"PEND", "STALE"}),
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_junit_xml(results: list[ConfigResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    total_tests = 0
    total_failures = 0
    total_skipped = 0
    testsuite = ET.Element("testsuite", name="vardict-rs-parity")

    for result in results:
        for chrom, chrom_result in result.chromosomes.items():
            total_tests += 1
            testcase = ET.SubElement(
                testsuite,
                "testcase",
                classname=f"parity.{result.config_id}",
                name=f"{result.label}.{chrom}",
                time=f"{result.elapsed_seconds:.3f}",
            )
            if chrom_result.status == "FAIL":
                total_failures += 1
                failure_message = "\n".join(_failure_lines(chrom_result.failures)) or "parity mismatch"
                ET.SubElement(testcase, "failure", message="parity mismatch").text = failure_message
            elif chrom_result.status in {"PEND", "STALE"}:
                total_skipped += 1
                ET.SubElement(testcase, "skipped", message=chrom_result.status.lower())

    testsuite.set("tests", str(total_tests))
    testsuite.set("failures", str(total_failures))
    testsuite.set("skipped", str(total_skipped))
    ET.indent(testsuite)
    ET.ElementTree(testsuite).write(path, encoding="utf-8", xml_declaration=True)


def _config_result_to_dict(result: ConfigResult) -> dict[str, object]:
    payload = asdict(result)
    payload["status"] = _config_status(result)
    return payload


def _config_status(result: ConfigResult) -> str:
    statuses = {chrom.status for chrom in result.chromosomes.values()}
    if "FAIL" in statuses:
        return "FAIL"
    if "STALE" in statuses:
        return "STALE"
    if "PEND" in statuses:
        return "PEND"
    if statuses == {"EMPTY"}:
        return "EMPTY"
    return "PASS"


def _failure_lines(failures: list[ShardFailure]) -> list[str]:
    lines: list[str] = []
    for failure in failures:
        detail = failure.summary or failure.reason
        if failure.first_diff_line is not None:
            detail = f"{detail} (line {failure.first_diff_line})"
        lines.append(f"{failure.shard_label} {failure.region}: {detail}")
    return lines


def _format_row(values: tuple[str, ...], widths: list[int]) -> str:
    return "  ".join(value.ljust(width) for value, width in zip(values, widths, strict=True))
