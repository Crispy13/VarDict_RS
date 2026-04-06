#!/usr/bin/env python3
"""Generate a self-contained HTML dashboard for parity sweep results."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import sys
from typing import Iterable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import parity_runner
from tests.parity.lib import config, reporter


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / parity_runner.BASE_DIR
DEFAULT_OUTPUT = PROJECT_ROOT / "tmp" / "parity_dashboard.html"

STATUS_PASS = "PASS"
STATUS_FAIL = "FAIL"
STATUS_EMPTY = "EMPTY"
STATUS_STALE = "STALE"
STATUS_PARTIAL = "PARTIAL"


@dataclass(frozen=True, slots=True)
class ConfigMeta:
    config_id: str
    label: str
    tier: int | str
    flags: str


@dataclass(frozen=True, slots=True)
class CellView:
    chrom: str
    result: reporter.ChromResult
    display_status: str
    tooltip: str


@dataclass(frozen=True, slots=True)
class SweepReportSummary:
    total_configs: int
    pass_configs: int
    fail_configs: int
    timestamp: str


@dataclass(frozen=True, slots=True)
class DashboardRow:
    meta: ConfigMeta
    result: reporter.ConfigResult
    cells: dict[str, CellView]


@dataclass(frozen=True, slots=True)
class DashboardModel:
    generated_at: str
    results_dir: Path
    rows: list[DashboardRow]
    chromosomes: list[str]
    sweep_report: SweepReportSummary | None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate a self-contained HTML parity dashboard.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Parity results directory (default: tmp/na12878_parity)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output HTML path (default: tmp/parity_dashboard.html)",
    )
    return parser


def fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def normalize_path(path: Path) -> Path:
    return path if path.is_absolute() else PROJECT_ROOT / path


def tier_sort_key(tier: int | str) -> tuple[int, str]:
    if isinstance(tier, int):
        return (0, f"{tier:03d}")
    return (1, str(tier))


def chromosome_sort_key(chrom: str) -> tuple[int, int | str]:
    normalized = chrom[3:] if chrom.startswith("chr") else chrom
    if normalized.isdigit():
        return (0, int(normalized))
    special = {"X": 23, "Y": 24, "MT": 25, "M": 25}
    if normalized in special:
        return (1, special[normalized])
    return (2, normalized)


def build_config_catalog() -> dict[str, ConfigMeta]:
    catalog: dict[str, ConfigMeta] = {}
    for entry in config.TEST_MATRIX:
        catalog[entry.label] = ConfigMeta(
            config_id=entry.config_id,
            label=entry.label,
            tier=entry.tier,
            flags=entry.flags,
        )

    pairwise_tsv = PROJECT_ROOT / config.DEFAULT_PAIRWISE_CONFIGS_TSV
    if pairwise_tsv.exists():
        for entry in config.load_pairwise_configs(pairwise_tsv):
            catalog.setdefault(
                entry.label,
                ConfigMeta(
                    config_id=entry.config_id,
                    label=entry.label,
                    tier=entry.tier,
                    flags=entry.flags,
                ),
            )
    return catalog


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_sweep_report(results_dir: Path) -> tuple[SweepReportSummary | None, dict[str, ConfigMeta], dict[str, dict[str, str]]]:
    report_path = results_dir / "option_parity_report.json"
    if not report_path.exists():
        return None, {}, {}

    payload = read_json(report_path)
    config_meta: dict[str, ConfigMeta] = {}
    config_statuses: dict[str, dict[str, str]] = {}

    for item in payload.get("configs", []):
        if not isinstance(item, dict):
            continue
        label = str(item.get("opts_label") or item.get("label") or item.get("config_id") or "unknown")
        config_meta[label] = ConfigMeta(
            config_id=str(item.get("config_id") or label),
            label=label,
            tier=item.get("tier", "?"),
            flags=str(item.get("cli_flags") or item.get("flags") or ""),
        )
        chromosomes = item.get("chromosomes", {})
        if isinstance(chromosomes, dict):
            config_statuses[label] = {str(chrom): str(status) for chrom, status in chromosomes.items()}

    summary = SweepReportSummary(
        total_configs=int(payload.get("total_configs", 0)),
        pass_configs=int(payload.get("pass_configs", 0)),
        fail_configs=int(payload.get("fail_configs", 0)),
        timestamp=str(payload.get("timestamp", "")),
    )
    return summary, config_meta, config_statuses


def parse_failures(payload: dict[str, object]) -> list[reporter.ShardFailure]:
    failures: list[reporter.ShardFailure] = []
    for item in payload.get("failures", []):
        if not isinstance(item, dict):
            continue
        failures.append(
            reporter.ShardFailure(
                shard_label=str(item.get("shard") or item.get("shard_label") or ""),
                region=str(item.get("region") or ""),
                reason=str(item.get("reason") or "unknown"),
                first_diff_line=_int_or_none(item.get("first_diff_line")),
                summary=str(item.get("summary") or ""),
            )
        )
    return failures


def parse_results_json(path: Path) -> tuple[reporter.ChromResult, float]:
    payload = read_json(path)
    result = reporter.ChromResult(
        status=str(payload.get("status", STATUS_PASS)),
        total_shards=int(payload.get("total_shards", 0)),
        pass_count=int(payload.get("pass", 0)),
        fail_count=int(payload.get("fail", 0)),
        empty_count=int(payload.get("empty", 0)),
        failures=parse_failures(payload),
    )
    return result, float(payload.get("elapsed_seconds", 0.0))


def _int_or_none(value: object) -> int | None:
    if value in {None, ""}:
        return None
    return int(value)


def load_marker_counts(marker_path: Path) -> tuple[int, int, int] | None:
    if not marker_path.exists():
        return None
    try:
        payload = read_json(marker_path)
    except json.JSONDecodeError:
        return None
    pass_count = int(payload.get("pass", 0))
    empty_count = int(payload.get("empty", 0))
    total_shards = int(payload.get("total_shards", pass_count + empty_count))
    return total_shards, pass_count, empty_count


def collect_status_markers(chrom_dir: Path) -> tuple[int, int, int, bool]:
    pass_count = 0
    fail_count = 0
    empty_count = 0
    saw_markers = False

    diff_dir = chrom_dir / "diff"
    if diff_dir.is_dir():
        for status_path in sorted(diff_dir.glob("shard_*.status")):
            saw_markers = True
            status = status_path.read_text(encoding="utf-8").strip().upper()
            if status == STATUS_FAIL:
                fail_count += 1
            elif status == STATUS_EMPTY:
                empty_count += 1
            elif status == STATUS_PASS:
                pass_count += 1
        for marker_path in sorted(diff_dir.glob("shard_*.verified")):
            saw_markers = True
            counts = load_marker_counts(marker_path)
            if counts is None:
                pass_count += 1
                continue
            _, marker_pass, marker_empty = counts
            pass_count += marker_pass or (0 if marker_empty else 1)
            empty_count += marker_empty

    fail_markers = list(chrom_dir.rglob("*.failed"))
    if fail_markers:
        saw_markers = True
        fail_count += len(fail_markers)

    has_artifacts = saw_markers or any((chrom_dir / name).exists() for name in ("java", "rust", "diff"))
    return pass_count, fail_count, empty_count, has_artifacts


def build_fallback_result(chrom_dir: Path, *, rust_mtime: str | None) -> tuple[reporter.ChromResult, str] | None:
    verified_path = chrom_dir / ".verified"
    if verified_path.exists():
        counts = load_marker_counts(verified_path) or (0, 0, 0)
        total_shards, pass_count, empty_count = counts
        display_status = STATUS_PASS if parity_runner.is_marker_current(verified_path, rust_mtime) else STATUS_STALE
        status = STATUS_EMPTY if total_shards > 0 and empty_count == total_shards else STATUS_PASS
        result = reporter.ChromResult(
            status=status,
            total_shards=total_shards,
            pass_count=pass_count,
            fail_count=0,
            empty_count=empty_count,
            failures=[],
        )
        return result, display_status

    pass_count, fail_count, empty_count, has_artifacts = collect_status_markers(chrom_dir)
    if not has_artifacts:
        return None

    total_shards = pass_count + fail_count + empty_count
    if fail_count > 0:
        status = STATUS_FAIL
        display_status = STATUS_FAIL
    elif total_shards == 0:
        status = STATUS_PARTIAL
        display_status = STATUS_PARTIAL
    elif empty_count == total_shards:
        status = STATUS_EMPTY
        display_status = STATUS_PARTIAL
    else:
        status = STATUS_PARTIAL
        display_status = STATUS_PARTIAL

    result = reporter.ChromResult(
        status=status,
        total_shards=total_shards,
        pass_count=pass_count,
        fail_count=fail_count,
        empty_count=empty_count,
        failures=[],
    )
    return result, display_status


def display_status_for_result(result: reporter.ChromResult, *, verified_path: Path, rust_mtime: str | None) -> str:
    if result.fail_count > 0 or result.status == STATUS_FAIL:
        return STATUS_FAIL
    if verified_path.exists() and not parity_runner.is_marker_current(verified_path, rust_mtime):
        return STATUS_STALE

    known_shards = result.pass_count + result.fail_count + result.empty_count
    if result.status == STATUS_EMPTY or (result.total_shards > 0 and result.empty_count == result.total_shards):
        return STATUS_EMPTY
    if result.status in {"PEND", STATUS_PARTIAL}:
        return STATUS_PARTIAL
    if result.total_shards > 0 and known_shards < result.total_shards:
        return STATUS_PARTIAL
    return STATUS_PASS


def failure_preview(failures: Iterable[reporter.ShardFailure], limit: int = 3) -> list[str]:
    lines: list[str] = []
    for failure in failures:
        detail = failure.summary or failure.reason
        if failure.first_diff_line is not None:
            detail = f"{detail} (line {failure.first_diff_line})"
        shard_label = failure.shard_label or "?"
        lines.append(f"{shard_label}: {detail}")
        if len(lines) >= limit:
            break
    return lines


def build_tooltip(meta: ConfigMeta, chrom: str, view: CellView) -> str:
    result = view.result
    lines = [
        f"{meta.config_id} [{meta.label}] chr{chrom}",
        f"Status: {view.display_status}",
        f"Shards: total={result.total_shards} pass={result.pass_count} fail={result.fail_count} empty={result.empty_count}",
    ]
    preview = failure_preview(result.failures)
    if preview:
        lines.append("Failures:")
        lines.extend(preview)
    return "\n".join(lines)


def collect_cell_view(chrom_dir: Path, chrom: str, meta: ConfigMeta, *, rust_mtime: str | None) -> CellView | None:
    results_path = chrom_dir / "results.json"
    if results_path.exists():
        result, _ = parse_results_json(results_path)
        display_status = display_status_for_result(result, verified_path=chrom_dir / ".verified", rust_mtime=rust_mtime)
        view = CellView(chrom=chrom, result=result, display_status=display_status, tooltip="")
        return CellView(chrom=chrom, result=result, display_status=display_status, tooltip=build_tooltip(meta, chrom, view))

    fallback = build_fallback_result(chrom_dir, rust_mtime=rust_mtime)
    if fallback is None:
        return None
    result, display_status = fallback
    view = CellView(chrom=chrom, result=result, display_status=display_status, tooltip="")
    return CellView(chrom=chrom, result=result, display_status=display_status, tooltip=build_tooltip(meta, chrom, view))


def fallback_meta(label: str, report_meta: dict[str, ConfigMeta], catalog: dict[str, ConfigMeta]) -> ConfigMeta:
    return report_meta.get(label) or catalog.get(label) or ConfigMeta(config_id=label, label=label, tier="?", flags="")


def collect_dashboard(results_dir: Path) -> DashboardModel:
    rust_mtime = parity_runner.rust_binary_mtime(parity_runner.DEFAULT_RUST_BIN)
    catalog = build_config_catalog()
    sweep_report, report_meta, report_statuses = parse_sweep_report(results_dir)

    labels: set[str] = set(report_meta)
    rows: list[DashboardRow] = []
    chromosomes: set[str] = set()

    if results_dir.exists():
        for child in results_dir.iterdir():
            if child.is_dir() and not child.name.startswith("."):
                labels.add(child.name)

    for label in sorted(labels, key=lambda value: (tier_sort_key(fallback_meta(value, report_meta, catalog).tier), fallback_meta(value, report_meta, catalog).config_id, value)):
        meta = fallback_meta(label, report_meta, catalog)
        label_dir = results_dir / label
        cell_views: dict[str, CellView] = {}
        elapsed_seconds = 0.0

        if label_dir.is_dir():
            for chrom_dir in sorted((path for path in label_dir.iterdir() if path.is_dir()), key=lambda path: chromosome_sort_key(path.name)):
                chrom = chrom_dir.name
                view = collect_cell_view(chrom_dir, chrom, meta, rust_mtime=rust_mtime)
                if view is None:
                    continue
                results_path = chrom_dir / "results.json"
                if results_path.exists():
                    _, cell_elapsed = parse_results_json(results_path)
                    elapsed_seconds += cell_elapsed
                cell_views[chrom] = view
                chromosomes.add(chrom)

        report_chrom_statuses = report_statuses.get(label, {})
        for chrom, status in sorted(report_chrom_statuses.items(), key=lambda item: chromosome_sort_key(item[0])):
            if chrom in cell_views:
                continue
            view = CellView(
                chrom=chrom,
                result=reporter.ChromResult(status=status, total_shards=0, pass_count=0, fail_count=0, empty_count=0, failures=[]),
                display_status=normalize_report_status(status),
                tooltip=f"{meta.config_id} [{meta.label}] chr{chrom}\nStatus: {normalize_report_status(status)}\nNo shard counts available in option_parity_report.json",
            )
            cell_views[chrom] = view
            chromosomes.add(chrom)

        if not cell_views:
            continue

        row_result = reporter.ConfigResult(
            config_id=meta.config_id,
            label=meta.label,
            tier=meta.tier,
            chromosomes={chrom: view.result for chrom, view in cell_views.items()},
            elapsed_seconds=elapsed_seconds,
        )
        rows.append(DashboardRow(meta=meta, result=row_result, cells=cell_views))

    rows.sort(key=lambda row: (tier_sort_key(row.meta.tier), row.meta.config_id, row.meta.label))

    return DashboardModel(
        generated_at=now_iso(),
        results_dir=results_dir,
        rows=rows,
        chromosomes=sorted(chromosomes, key=chromosome_sort_key),
        sweep_report=sweep_report,
    )


def normalize_report_status(status: str) -> str:
    normalized = status.upper()
    if normalized in {STATUS_PASS, STATUS_FAIL, STATUS_EMPTY, STATUS_STALE}:
        return normalized
    return STATUS_PARTIAL


def summarize_statuses(statuses: Iterable[str]) -> str:
    values = list(statuses)
    if not values:
        return STATUS_EMPTY
    unique = set(values)
    if STATUS_FAIL in unique:
        return STATUS_FAIL
    if STATUS_PARTIAL in unique:
        return STATUS_PARTIAL
    if unique == {STATUS_STALE}:
        return STATUS_STALE
    if STATUS_STALE in unique:
        return STATUS_STALE
    if unique == {STATUS_EMPTY}:
        return STATUS_EMPTY
    if unique <= {STATUS_PASS, STATUS_EMPTY}:
        return STATUS_PASS if STATUS_PASS in unique else STATUS_EMPTY
    return STATUS_PARTIAL


def status_css_class(status: str) -> str:
    return {
        STATUS_PASS: "pass",
        STATUS_FAIL: "fail",
        STATUS_EMPTY: "empty",
        STATUS_STALE: "stale",
        STATUS_PARTIAL: "partial",
    }.get(status, "partial")


def percent(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.0%"
    return f"{(100.0 * numerator / denominator):.1f}%"


def row_cell_counts(row: DashboardRow) -> dict[str, int]:
    counts = {STATUS_PASS: 0, STATUS_FAIL: 0, STATUS_EMPTY: 0, STATUS_STALE: 0, STATUS_PARTIAL: 0}
    for view in row.cells.values():
        counts[view.display_status] += 1
    return counts


def chrom_cell_counts(model: DashboardModel, chrom: str) -> dict[str, int]:
    counts = {STATUS_PASS: 0, STATUS_FAIL: 0, STATUS_EMPTY: 0, STATUS_STALE: 0, STATUS_PARTIAL: 0}
    for row in model.rows:
        view = row.cells.get(chrom)
        if view is None:
            continue
        counts[view.display_status] += 1
    return counts


def total_shard_counts(model: DashboardModel) -> tuple[int, int, int, int]:
    total = 0
    passing = 0
    failing = 0
    empty = 0
    for row in model.rows:
        for view in row.cells.values():
            total += view.result.total_shards
            passing += view.result.pass_count
            failing += view.result.fail_count
            empty += view.result.empty_count
    return total, passing, failing, empty


def overall_cell_counts(model: DashboardModel) -> dict[str, int]:
    counts = {STATUS_PASS: 0, STATUS_FAIL: 0, STATUS_EMPTY: 0, STATUS_STALE: 0, STATUS_PARTIAL: 0}
    for row in model.rows:
        row_counts = row_cell_counts(row)
        for status, value in row_counts.items():
            counts[status] += value
    return counts


def legend_html() -> str:
    items = [
        (STATUS_PASS, "PASS"),
        (STATUS_FAIL, "FAIL"),
        (STATUS_PARTIAL, "PARTIAL"),
        (STATUS_EMPTY, "EMPTY"),
        (STATUS_STALE, "STALE"),
    ]
    parts = []
    for status, label in items:
        css_class = status_css_class(status)
        parts.append(
            f'<div class="legend-item"><span class="legend-swatch {css_class}"></span><span>{html.escape(label)}</span></div>'
        )
    return "".join(parts)


def render_summary_cards(model: DashboardModel) -> str:
    total, passing, failing, empty = total_shard_counts(model)
    cell_counts = overall_cell_counts(model)
    banner_status = summarize_statuses(
        status
        for status, count in cell_counts.items()
        for _ in range(count)
    )

    cards = [
        (
            "Overall shard parity",
            f"{passing}/{total} shards passing ({percent(passing, total)})" if total else "0/0 shards passing (0.0%)",
            f"{failing} failing, {empty} empty",
        ),
        (
            "Cells",
            f"{sum(cell_counts.values())} populated cells",
            f"{cell_counts[STATUS_PASS]} pass, {cell_counts[STATUS_FAIL]} fail, {cell_counts[STATUS_PARTIAL]} partial",
        ),
        (
            "Configs",
            f"{len(model.rows)} configs",
            f"{len(model.chromosomes)} chromosomes",
        ),
    ]

    if model.sweep_report is not None:
        cards.append(
            (
                "Sweep report",
                f"{model.sweep_report.pass_configs}/{model.sweep_report.total_configs} configs passing",
                f"{model.sweep_report.fail_configs} failing, generated {model.sweep_report.timestamp or 'unknown'}",
            )
        )

    pieces = [
        f'<section class="banner {status_css_class(banner_status)}"><div class="banner-title">VarDict-rs Parity Dashboard</div><div class="banner-subtitle">Generated {html.escape(model.generated_at)} from {html.escape(model.results_dir.as_posix())}</div></section>',
        '<section class="cards">',
    ]
    for title, primary, secondary in cards:
        pieces.append(
            '<article class="card">'
            f'<div class="card-title">{html.escape(title)}</div>'
            f'<div class="card-primary">{html.escape(primary)}</div>'
            f'<div class="card-secondary">{html.escape(secondary)}</div>'
            '</article>'
        )
    pieces.append("</section>")
    return "".join(pieces)


def render_grid(model: DashboardModel) -> str:
    if not model.rows or not model.chromosomes:
        return (
            '<section class="empty-state">'
            '<h2>No results found</h2>'
            '<p>No parity result cells were discovered in the selected results directory.</p>'
            '</section>'
        )

    header_cells = ['<th class="sticky-col header-col">Config</th>']
    for chrom in model.chromosomes:
        header_cells.append(f'<th class="header-chrom">chr{html.escape(chrom)}</th>')
    header_cells.append('<th class="summary-col">Summary</th>')

    body_rows: list[str] = []
    for row in model.rows:
        counts = row_cell_counts(row)
        row_status = summarize_statuses(view.display_status for view in row.cells.values())
        row_summary = (
            f"{counts[STATUS_PASS]} pass, {counts[STATUS_FAIL]} fail, "
            f"{counts[STATUS_PARTIAL]} partial, {counts[STATUS_EMPTY] + counts[STATUS_STALE]} gray"
        )
        row_title = html.escape(
            f"{row.meta.config_id} | tier {row.meta.tier} | {row.meta.flags or 'no extra flags'}"
        )
        cells = [
            '<th class="sticky-col row-label" '
            f'title="{row_title}">'
            f'<span class="config-id">{html.escape(row.meta.config_id)}</span>'
            f'<span class="config-label">{html.escape(row.meta.label)}</span>'
            '</th>'
        ]
        for chrom in model.chromosomes:
            view = row.cells.get(chrom)
            if view is None:
                cells.append('<td class="matrix-cell missing"><span class="swatch ghost"></span></td>')
                continue
            css_class = status_css_class(view.display_status)
            aria = html.escape(f"{row.meta.label} chr{chrom}: {view.display_status}")
            tooltip = html.escape(view.tooltip)
            cells.append(
                f'<td class="matrix-cell {css_class}" title="{tooltip}">' 
                f'<span class="swatch {css_class}" aria-label="{aria}"></span>'
                '</td>'
            )

        cells.append(
            f'<td class="summary {status_css_class(row_status)}" title="{html.escape(row_summary)}">'
            f'<div class="summary-status">{html.escape(row_status)}</div>'
            f'<div class="summary-text">{html.escape(row_summary)}</div>'
            '</td>'
        )
        body_rows.append(f"<tr>{''.join(cells)}</tr>")

    summary_cells = ['<th class="sticky-col footer-label">Chrom summary</th>']
    for chrom in model.chromosomes:
        counts = chrom_cell_counts(model, chrom)
        chrom_status = summarize_statuses(
            row.cells[chrom].display_status
            for row in model.rows
            if chrom in row.cells
        )
        chrom_summary = (
            f"{counts[STATUS_PASS]} pass, {counts[STATUS_FAIL]} fail, "
            f"{counts[STATUS_PARTIAL]} partial, {counts[STATUS_EMPTY] + counts[STATUS_STALE]} gray"
        )
        summary_cells.append(
            f'<td class="summary {status_css_class(chrom_status)}" title="{html.escape(chrom_summary)}">'
            f'<div class="summary-status">{html.escape(chrom_status)}</div>'
            f'<div class="summary-text">{html.escape(chrom_summary)}</div>'
            '</td>'
        )
    summary_cells.append('<td class="summary footer-corner">Config totals</td>')

    return (
        '<section class="legend">'
        f'{legend_html()}'
        '</section>'
        '<section class="matrix-wrap">'
        '<table class="matrix">'
        f'<thead><tr>{"".join(header_cells)}</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody>'
        f'<tfoot><tr>{"".join(summary_cells)}</tr></tfoot>'
        '</table>'
        '</section>'
    )


def render_toolbar() -> str:
    return (
        '<div class="refresh-overlay" id="refresh-overlay">'
        '<div class="refresh-overlay-box">'
        '<div class="refresh-spinner"></div>'
        '<div class="refresh-msg" id="refresh-msg">Refreshing dashboard…</div>'
        '</div>'
        '</div>'
        '<section class="toolbar">'
        '<button class="btn-refresh" id="btn-refresh" onclick="doRefresh()">&#8635; Refresh</button>'
        '<span class="auto-label">Auto-refresh:</span>'
        '<button class="btn-auto active" id="auto-off" onclick="setAuto(0)">Off</button>'
        '<button class="btn-auto" id="auto-30" onclick="setAuto(30)">30 s</button>'
        '<button class="btn-auto" id="auto-60" onclick="setAuto(60)">60 s</button>'
        '<span class="countdown" id="countdown"></span>'
        '</section>'
    )


REFRESH_JS = """
  <script>
    var _tickTimer = null;
    var _remaining = 0;

    function showOverlay(msg) {
      document.getElementById('refresh-msg').textContent = msg || 'Refreshing dashboard\u2026';
      document.getElementById('refresh-overlay').style.display = 'flex';
    }

    function doRefresh() {
      showOverlay('Refreshing dashboard\u2026');
      // Small delay so the overlay renders before the browser starts loading
      setTimeout(function() { location.reload(); }, 60);
    }

    function setAuto(seconds) {
      clearInterval(_tickTimer);
      document.querySelectorAll('.btn-auto').forEach(function(b) { b.classList.remove('active'); });
      var id = seconds === 0 ? 'auto-off' : 'auto-' + seconds;
      var el = document.getElementById(id);
      if (el) el.classList.add('active');
      document.getElementById('countdown').textContent = '';
      if (seconds > 0) {
        _remaining = seconds;
        document.getElementById('countdown').textContent = 'next in ' + _remaining + 's';
        _tickTimer = setInterval(function() {
          _remaining -= 1;
          if (_remaining <= 0) {
            clearInterval(_tickTimer);
            showOverlay('Auto-refreshing\u2026');
            setTimeout(function() { location.reload(); }, 60);
          } else {
            document.getElementById('countdown').textContent = 'next in ' + _remaining + 's';
          }
        }, 1000);
      }
    }

    // Show overlay while page is loading (covers the flash of old content)
    window.addEventListener('beforeunload', function() {
      showOverlay('Loading new data\u2026');
    });
  </script>
"""


def render_html(model: DashboardModel) -> str:
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>VarDict-rs Parity Dashboard</title>
  <style>
    :root {{
      --bg: #f6f3ec;
      --panel: #fffdf8;
      --ink: #1c1a17;
      --muted: #6d655a;
      --grid: #d8cfbf;
      --pass: #3f8f5c;
      --fail: #bf3f3f;
      --partial: #d1a52a;
      --empty: #9a9488;
      --stale: #7e7a73;
      --shadow: rgba(28, 26, 23, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(209, 165, 42, 0.12), transparent 30%),
        linear-gradient(180deg, #fbf7ef 0%, var(--bg) 100%);
      font-family: "Iowan Old Style", "Palatino Linotype", "Book Antiqua", Georgia, serif;
    }}
    .page {{
      padding: 24px;
      display: grid;
      gap: 18px;
    }}
    .banner {{
      padding: 20px 22px;
      border: 1px solid var(--grid);
      border-left-width: 8px;
      background: var(--panel);
      box-shadow: 0 14px 30px var(--shadow);
    }}
    .banner.pass {{ border-left-color: var(--pass); }}
    .banner.fail {{ border-left-color: var(--fail); }}
    .banner.partial {{ border-left-color: var(--partial); }}
    .banner.empty, .banner.stale {{ border-left-color: var(--empty); }}
    .banner-title {{ font-size: 1.8rem; font-weight: 700; letter-spacing: 0.02em; }}
    .banner-subtitle {{ margin-top: 6px; color: var(--muted); font-size: 0.98rem; }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 14px;
    }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--grid);
      box-shadow: 0 10px 24px var(--shadow);
      padding: 16px 18px;
    }}
    .card-title {{ font-size: 0.8rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.08em; }}
    .card-primary {{ margin-top: 8px; font-size: 1.35rem; font-weight: 700; }}
    .card-secondary {{ margin-top: 6px; color: var(--muted); font-size: 0.95rem; }}
    .legend {{
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      align-items: center;
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      font-size: 0.92rem;
    }}
    .legend-item {{ display: inline-flex; align-items: center; gap: 8px; }}
    .legend-swatch, .swatch {{
      width: 16px;
      height: 16px;
      border: 1px solid rgba(0, 0, 0, 0.2);
      display: inline-block;
      border-radius: 3px;
    }}
    .pass {{ background: var(--pass); }}
    .fail {{ background: var(--fail); }}
    .partial {{ background: var(--partial); }}
    .empty {{ background: var(--empty); }}
    .stale {{ background: var(--stale); }}
    .ghost {{ background: repeating-linear-gradient(135deg, #f4eee2, #f4eee2 4px, #ece4d8 4px, #ece4d8 8px); }}
    .matrix-wrap {{
      overflow: auto;
      border: 1px solid var(--grid);
      background: var(--panel);
      box-shadow: 0 12px 28px var(--shadow);
      max-height: calc(100vh - 260px);
    }}
    .matrix {{
      border-collapse: separate;
      border-spacing: 0;
      width: max-content;
      min-width: 100%;
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      font-size: 0.85rem;
    }}
    .matrix th, .matrix td {{
      border-right: 1px solid var(--grid);
      border-bottom: 1px solid var(--grid);
      padding: 8px;
      background: var(--panel);
      text-align: center;
      vertical-align: middle;
    }}
    .matrix thead th {{
      position: sticky;
      top: 0;
      z-index: 3;
      background: #f4ecdd;
    }}
    .sticky-col {{
      position: sticky;
      left: 0;
      z-index: 2;
      background: #fbf6ea;
      text-align: left;
      min-width: 180px;
      max-width: 180px;
    }}
    .matrix thead .sticky-col {{ z-index: 4; }}
    .header-col, .footer-label {{ text-transform: uppercase; letter-spacing: 0.06em; font-size: 0.78rem; }}
    .header-chrom {{ min-width: 48px; }}
    .row-label {{ line-height: 1.25; }}
    .config-id {{ display: block; font-weight: 700; }}
    .config-label {{ display: block; color: var(--muted); font-size: 0.82rem; }}
    .matrix-cell {{ width: 40px; min-width: 40px; padding: 6px; }}
    .matrix-cell .swatch {{ width: 20px; height: 20px; }}
    .summary-col {{ min-width: 180px; }}
    .summary {{ min-width: 180px; text-align: left; }}
    .summary-status {{ font-weight: 700; }}
    .summary-text {{ margin-top: 4px; color: var(--muted); font-size: 0.78rem; line-height: 1.35; }}
    .footer-corner {{ color: var(--muted); font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.06em; }}
    .empty-state {{
      background: var(--panel);
      border: 1px dashed var(--grid);
      box-shadow: 0 10px 24px var(--shadow);
      padding: 28px;
      text-align: center;
    }}
    .empty-state h2 {{ margin: 0 0 10px; }}
    .empty-state p {{ margin: 0; color: var(--muted); }}
    @media (max-width: 900px) {{
      .page {{ padding: 14px; }}
      .banner-title {{ font-size: 1.45rem; }}
      .sticky-col {{ min-width: 150px; max-width: 150px; }}
      .summary, .summary-col {{ min-width: 150px; }}
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      gap: 10px;
      flex-wrap: wrap;
      background: var(--panel);
      border: 1px solid var(--grid);
      padding: 10px 16px;
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      font-size: 0.88rem;
      box-shadow: 0 6px 16px var(--shadow);
    }}
    .btn-refresh {{
      padding: 6px 14px;
      background: var(--ink);
      color: #fffdf8;
      border: none;
      border-radius: 4px;
      cursor: pointer;
      font-size: 0.9rem;
      font-family: inherit;
    }}
    .btn-refresh:hover {{ opacity: 0.82; }}
    .auto-label {{ color: var(--muted); }}
    .btn-auto {{
      padding: 5px 11px;
      background: var(--panel);
      color: var(--ink);
      border: 1px solid var(--grid);
      border-radius: 4px;
      cursor: pointer;
      font-size: 0.85rem;
      font-family: inherit;
    }}
    .btn-auto:hover {{ border-color: var(--ink); }}
    .btn-auto.active {{
      background: var(--ink);
      color: #fffdf8;
      border-color: var(--ink);
    }}
    .countdown {{ color: var(--muted); min-width: 90px; }}
    .refresh-overlay {{
      display: none;
      position: fixed;
      inset: 0;
      z-index: 9999;
      background: rgba(246, 243, 236, 0.82);
      backdrop-filter: blur(3px);
      align-items: center;
      justify-content: center;
    }}
    .refresh-overlay-box {{
      background: var(--panel);
      border: 1px solid var(--grid);
      box-shadow: 0 20px 48px var(--shadow);
      padding: 36px 52px;
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 18px;
    }}
    .refresh-spinner {{
      width: 36px;
      height: 36px;
      border: 4px solid var(--grid);
      border-top-color: var(--ink);
      border-radius: 50%;
      animation: spin 0.7s linear infinite;
    }}
    @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
    .refresh-msg {{
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      font-size: 1rem;
      color: var(--ink);
    }}
    @media print {{
      body {{ background: white; }}
      .page {{ padding: 0; gap: 10px; }}
      .banner, .card, .matrix-wrap, .empty-state {{ box-shadow: none; }}
      .matrix-wrap {{ max-height: none; overflow: visible; }}
      .toolbar {{ display: none; }}
    }}
  </style>
</head>
<body>
  <main class="page">
    {render_toolbar()}
    {render_summary_cards(model)}
    {render_grid(model)}
  </main>
  {REFRESH_JS}
</body>
</html>
'''


def write_dashboard(output_path: Path, model: DashboardModel) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_html(model), encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    results_dir = normalize_path(args.results_dir)
    output_path = normalize_path(args.output)

    if results_dir.exists() and not results_dir.is_dir():
        return fail(f"Results path is not a directory: {results_dir}")

    model = collect_dashboard(results_dir)
    write_dashboard(output_path, model)

    print(f"Wrote parity dashboard to {output_path}")
    if not model.rows or not model.chromosomes:
        print("No results found")
    else:
        total, passing, failing, empty = total_shard_counts(model)
        print(
            f"Included {len(model.rows)} configs x {len(model.chromosomes)} chromosomes; "
            f"shards total={total} pass={passing} fail={failing} empty={empty}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())