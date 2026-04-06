"""Shard generation helpers for the parity harness."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class Shard:
    index: int
    label: str
    chrom: str
    start: int
    end: int


def load_fai(fai_path: Path) -> dict[str, int]:
    chromosomes: dict[str, int] = {}
    with fai_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chromosomes[fields[0]] = int(fields[1])
    return chromosomes


def generate_shards(chrom: str, chrom_len: int, shard_size: int = 1_000_000) -> list[Shard]:
    if chrom_len <= 0:
        return []
    if shard_size <= 0:
        raise ValueError("shard_size must be greater than zero")

    shards: list[Shard] = []
    index = 0
    start = 1
    while start <= chrom_len:
        end = min(start + shard_size - 1, chrom_len)
        shards.append(
            Shard(
                index=index,
                label=f"{index + 1:03d}",
                chrom=chrom,
                start=start,
                end=end,
            )
        )
        index += 1
        start = end + 1
    return shards


def shard_region_string(shard: Shard) -> str:
    return f"{shard.chrom}:{shard.start}-{shard.end}"


def shard_bed_line(shard: Shard) -> str:
    return f"{shard.chrom}\t{shard.start - 1}\t{shard.end}"
