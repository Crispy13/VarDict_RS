"""Parity harness config definitions mirrored from the bash scripts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import shlex

Tier = int | str

DEFAULT_PAIRWISE_CONFIGS_TSV = Path("tests/pairwise_configs.tsv")


@dataclass(frozen=True, slots=True)
class TestConfig:
    config_id: str
    label: str
    flags: str
    tier: Tier
    blocked: bool = False

    @property
    def uses_pileup(self) -> bool:
        return "-p" in shlex.split(self.flags)


@dataclass(frozen=True, slots=True)
class Preset:
    name: str
    config_ids: tuple[str, ...] = ()
    tier_filters: tuple[Tier, ...] = ()
    chromosomes: tuple[str, ...] = ()
    all_real_chromosomes: bool = False
    include_pairwise: bool = False
    stages: tuple[str, ...] = ()


TEST_MATRIX: list[TestConfig] = [
    TestConfig("T1-03", "nosv", "-U", 1),
    TestConfig("T1-05", "freq-low", "-f 0.001", 1),
    TestConfig("T1-13", "fisher", "--fisher", 1),
    TestConfig("T1-04", "nosv-dedup", "-U --deldupvar", 1),
    TestConfig("T1-06", "freq-high", "-f 0.05", 1),
    TestConfig("T1-07", "freq-zero", "-f 0.0", 1),
    TestConfig("T1-08", "minr-1", "-r 1", 1),
    TestConfig("T1-09", "minr-5", "-r 5", 1),
    TestConfig("T1-11", "mapq-10", "-Q 10", 1),
    TestConfig("T1-12", "mapq-30", "-Q 30", 1),
    TestConfig("T2-09", "filter-0x700", "-F 0x700", 2),
    TestConfig("T2-16", "chimeric", "--chimeric", 2),
    TestConfig("T2-01", "qual-0", "-q 0", 2),
    TestConfig("T2-02", "qual-10", "-q 10", 2),
    TestConfig("T2-03", "qual-30", "-q 30", 2),
    TestConfig("T2-04", "qual-40", "-q 40", 2),
    TestConfig("T2-05", "vext-0", "-X 0", 2),
    TestConfig("T2-06", "vext-5", "-X 5", 2),
    TestConfig("T2-07", "mismatch-3", "-m 3", 2),
    TestConfig("T2-08", "mismatch-15", "-m 15", 2),
    TestConfig("T2-10", "filter-0x100", "-F 0x100", 2),
    TestConfig("T2-11", "indel-3prime", "-3", 2),
    TestConfig("T2-12", "readpos-0", "-P 0", 2),
    TestConfig("T2-13", "readpos-10", "-P 10", 2),
    TestConfig("T2-14", "qratio-low", "-o 0.5", 2),
    TestConfig("T2-15", "qratio-high", "-o 2.0", 2),
    TestConfig("T3-01", "clinical-wgs", "-f 0.001 -Q 10 -F 0x700", 3),
    TestConfig("T4-04", "inssize-small", "-w 200 -W 50", 4),
    TestConfig("T3-04", "fisher-lowfreq", "-f 0.001 --fisher", 3),
    TestConfig("T3-05", "nosv-tight", "-U -f 0.001 -Q 10", 3),
    TestConfig("T3-06", "3prime-lowfreq", "-3 -f 0.001", 3),
    TestConfig("T3-08", "quality-gauntlet", "-M 25 -m 5 -Q 10", 3),
    TestConfig("T4-01", "extend-150", "-x 150", 4),
    TestConfig("T4-02", "refext-600", "-Y 600", 4),
    TestConfig("T4-03", "refext-2000", "-Y 2000", 4),
    TestConfig("T4-05", "trim-130", "-T 130", 4),
    TestConfig("T1-10", "no-realign", "-k 0", 1),
    TestConfig("T1-14", "debug", "-D", 1),
    TestConfig("T3-03", "no-realign-lowfreq", "-k 0 -f 0.001", 3),
    TestConfig("T1-01", "pileup", "-p", 1),
    TestConfig("T1-02", "pileup-max", "-p -f 0.0 -r 1", 1),
    TestConfig("T3-02", "pileup-strict", "-p -q 30", 3),
    TestConfig("T3-07", "pileup-relaxed", "-p -r 1 -q 10 -Q 0", 3),
    TestConfig("T4-06", "include-n-pileup", "-K -p -r 1 -f 0.0", 4),
]

TEST_MATRIX_BY_ID: dict[str, TestConfig] = {config.config_id: config for config in TEST_MATRIX}

PRESETS: dict[str, Preset] = {
    "smoke": Preset(
        name="smoke",
        config_ids=("T1-01", "T1-03", "T1-13"),
        chromosomes=("20", "22", "MT"),
    ),
    "dev": Preset(
        name="dev",
        config_ids=("T1-01", "T1-03", "T1-05", "T1-10", "T1-13", "T1-14", "T2-09", "T2-16", "T3-01", "T4-04"),
        chromosomes=("1", "2", "5", "10", "14", "17", "20", "22", "X", "MT"),
    ),
    "tier1": Preset(
        name="tier1",
        tier_filters=(1,),
        chromosomes=("20", "22", "MT"),
    ),
    "config-spread": Preset(
        name="config-spread",
        chromosomes=("20", "22", "MT"),
    ),
    "core-wide": Preset(
        name="core-wide",
        tier_filters=(1,),
        all_real_chromosomes=True,
    ),
    "pairwise": Preset(
        name="pairwise",
        tier_filters=("PW",),
        chromosomes=("1", "2", "5", "10", "14", "17", "20", "22", "X", "MT"),
        include_pairwise=True,
    ),
    "full-gate": Preset(
        name="full-gate",
        stages=("smoke", "tier1", "config-spread", "core-wide"),
    ),
    "release": Preset(
        name="release",
        all_real_chromosomes=True,
    ),
}


def load_pairwise_configs(tsv_path: Path = DEFAULT_PAIRWISE_CONFIGS_TSV) -> list[TestConfig]:
    pairwise_configs: list[TestConfig] = []
    with tsv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            label = (row.get("label") or "").strip()
            if not label:
                continue
            pairwise_configs.append(
                TestConfig(
                    config_id=label,
                    label=label,
                    flags=(row.get("cli_flags") or "").strip(),
                    tier="PW",
                )
            )
    return pairwise_configs


def resolve_preset(
    name: str,
    *,
    pairwise_tsv_path: Path = DEFAULT_PAIRWISE_CONFIGS_TSV,
    real_chromosomes: list[str] | tuple[str, ...] | None = None,
) -> tuple[list[TestConfig], list[str]]:
    preset = PRESETS[name]
    if preset.stages:
        raise ValueError(f"Preset '{name}' is sequential; use resolve_preset_sequence().")

    configs: list[TestConfig]
    if preset.include_pairwise:
        configs = load_pairwise_configs(pairwise_tsv_path)
    elif preset.config_ids:
        configs = [TEST_MATRIX_BY_ID[config_id] for config_id in preset.config_ids]
    elif preset.tier_filters:
        allowed_tiers = set(preset.tier_filters)
        configs = [config for config in TEST_MATRIX if config.tier in allowed_tiers]
    else:
        configs = list(TEST_MATRIX)

    if preset.all_real_chromosomes:
        if real_chromosomes is None:
            raise ValueError(f"Preset '{name}' requires the real chromosome list from the FAI.")
        chromosomes = list(real_chromosomes)
    else:
        chromosomes = list(preset.chromosomes)

    return configs, chromosomes


def resolve_preset_sequence(
    name: str,
    *,
    pairwise_tsv_path: Path = DEFAULT_PAIRWISE_CONFIGS_TSV,
    real_chromosomes: list[str] | tuple[str, ...] | None = None,
) -> list[tuple[list[TestConfig], list[str]]]:
    preset = PRESETS[name]
    if not preset.stages:
        return [resolve_preset(name, pairwise_tsv_path=pairwise_tsv_path, real_chromosomes=real_chromosomes)]
    return [
        resolve_preset(stage_name, pairwise_tsv_path=pairwise_tsv_path, real_chromosomes=real_chromosomes)
        for stage_name in preset.stages
    ]


def resource_overrides(config: TestConfig) -> dict[str, int | str]:
    if config.uses_pileup:
        return {
            "java_parallel": 5,
            "java_heap": "4g",
            "parallel": 10,
        }
    return {}
