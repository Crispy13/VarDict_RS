"""Resource monitoring helpers for the parity harness."""

from __future__ import annotations

from pathlib import Path
import shutil


def mem_available_gb() -> float:
    with Path("/proc/meminfo").open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                kb = int(line.split()[1])
                return kb / 1024 / 1024
    raise RuntimeError("MemAvailable not found in /proc/meminfo")


def disk_available_gb(path: Path = Path("tmp")) -> float:
    usage = shutil.disk_usage(path)
    return usage.free / 1024 / 1024 / 1024


def check_resources(
    warn_mem_gb: float,
    abort_mem_gb: float,
    warn_disk_gb: float,
    *,
    path: Path = Path("tmp"),
    abort_disk_gb: float = 5.0,
) -> list[str]:
    messages: list[str] = []
    mem_gb = mem_available_gb()
    disk_gb = disk_available_gb(path)

    if mem_gb < abort_mem_gb:
        messages.append(f"ABORT: MemAvailable {mem_gb:.1f}GB < {abort_mem_gb:.1f}GB")
    elif mem_gb < warn_mem_gb:
        messages.append(f"WARN: MemAvailable {mem_gb:.1f}GB < {warn_mem_gb:.1f}GB")

    if disk_gb < abort_disk_gb:
        messages.append(f"ABORT: disk {disk_gb:.1f}GB < {abort_disk_gb:.1f}GB")
    elif disk_gb < warn_disk_gb:
        messages.append(f"WARN: disk {disk_gb:.1f}GB < {warn_disk_gb:.1f}GB")

    return messages