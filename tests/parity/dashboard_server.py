#!/usr/bin/env python3
"""Live HTTP server for the parity dashboard.

Serves the dashboard HTML dynamically — every page load re-reads ``results_dir``
and regenerates fresh content without writing any intermediate file.

Usage::

    python -m tests.parity serve [--results-dir tmp/na12878_parity] [--port 7777]
"""

from __future__ import annotations

import argparse
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
from pathlib import Path
import sys
import threading
import webbrowser
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import dashboard as _dashboard
from tests.parity import parity_runner as _runner


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / _runner.BASE_DIR
DEFAULT_PORT = 7777


def _build_handler_class(results_dir: Path) -> type[BaseHTTPRequestHandler]:
    class DashboardHandler(BaseHTTPRequestHandler):
        _results_dir: Path = results_dir

        def log_message(self, fmt: str, *args: object) -> None:  # noqa: ANN001
            # Suppress default per-request noise; errors still go to log_error
            pass

        def do_GET(self) -> None:  # noqa: N802
            path = self.path.split("?")[0]
            if path in {"/", "/index.html"}:
                self._serve_dashboard()
            elif path == "/ping":
                self._serve_json({"status": "ok"})
            else:
                self.send_error(404, "Not found")

        def _serve_dashboard(self) -> None:
            model = _dashboard.collect_dashboard(self._results_dir)
            body = _dashboard.render_html(model).encode("utf-8")
            total, passing, failing, empty = _dashboard.total_shard_counts(model)
            print(
                f"[dashboard] {len(model.rows)} configs × {len(model.chromosomes)} chromosomes | "
                f"shards total={total} pass={passing} fail={failing} empty={empty}"
            )
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _serve_json(self, payload: dict[str, object]) -> None:
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

    return DashboardHandler


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Serve the parity dashboard over HTTP (live).")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help=f"Parity results directory (default: {DEFAULT_RESULTS_DIR.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to listen on (default: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--no-open",
        action="store_true",
        help="Don't open browser automatically",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    results_dir = args.results_dir if args.results_dir.is_absolute() else PROJECT_ROOT / args.results_dir

    handler_class = _build_handler_class(results_dir)
    url = f"http://localhost:{args.port}/"

    try:
        server = HTTPServer(("localhost", args.port), handler_class)
    except OSError as exc:
        print(f"ERROR: Cannot bind port {args.port}: {exc}", file=sys.stderr)
        return 1

    print(f"Serving parity dashboard at  {url}")
    print(f"Results directory            {results_dir}")
    print(f"Refresh button regenerates   yes (live, no caching)")
    print("Press Ctrl+C to stop.\n")

    if not args.no_open:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
