#!/usr/bin/env python3
"""Unified entry point for the parity harness tools."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Callable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import dashboard, dashboard_server, golden_archive, option_parity_runner, parity_cleanup, parity_runner, parity_status


CommandHandler = Callable[[Sequence[str] | None], int]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Unified CLI for VarDict parity tools.")
    subparsers = parser.add_subparsers(dest="command")
    subparsers.add_parser("run", add_help=False, help="Run single-config parity via parity_runner.py")
    subparsers.add_parser("sweep", add_help=False, help="Run config sweep via option_parity_runner.py")
    subparsers.add_parser("status", add_help=False, help="Display parity status matrix")
    subparsers.add_parser("cleanup", add_help=False, help="Clean stale or scoped parity artifacts")
    subparsers.add_parser("archive", add_help=False, help="Create or restore Java golden archives")
    subparsers.add_parser("dashboard", add_help=False, help="Generate a static HTML parity dashboard")
    subparsers.add_parser("serve", add_help=False, help="Serve live parity dashboard over HTTP (with working Refresh)")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args, remaining = parser.parse_known_args(argv)
    if args.command is None:
        parser.print_help()
        return 0

    handlers: dict[str, CommandHandler] = {
        "run": parity_runner.main,
        "sweep": option_parity_runner.main,
        "status": parity_status.main,
        "cleanup": parity_cleanup.main,
        "archive": golden_archive.main,
        "dashboard": dashboard.main,
        "serve": dashboard_server.main,
    }
    return handlers[args.command](remaining)


if __name__ == "__main__":
    raise SystemExit(main())