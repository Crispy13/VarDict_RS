"""Checkpoint helpers for parity sweep resume support."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
import json

CHECKPOINT_PATH = Path("tmp/na12878_parity/.sweep_checkpoint.json")


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(slots=True)
class SweepCheckpoint:
    completed: list[tuple[str, str]] = field(default_factory=list)
    failed: list[tuple[str, str, str]] = field(default_factory=list)
    started_at: str = field(default_factory=_timestamp)
    last_updated: str = field(default_factory=_timestamp)


def load_checkpoint(path: Path = CHECKPOINT_PATH) -> SweepCheckpoint | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    return SweepCheckpoint(
        completed=[tuple(item) for item in payload.get("completed", [])],
        failed=[tuple(item) for item in payload.get("failed", [])],
        started_at=payload.get("started_at", _timestamp()),
        last_updated=payload.get("last_updated", _timestamp()),
    )


def save_checkpoint(checkpoint: SweepCheckpoint, path: Path = CHECKPOINT_PATH) -> None:
    checkpoint.last_updated = _timestamp()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "completed": checkpoint.completed,
        "failed": checkpoint.failed,
        "started_at": checkpoint.started_at,
        "last_updated": checkpoint.last_updated,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def clear_checkpoint(path: Path = CHECKPOINT_PATH) -> None:
    path.unlink(missing_ok=True)
