"""Java and Rust shard runners for the parity harness."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import gzip
import shutil
import shlex
import subprocess
import tempfile
import time

from .config import TestConfig
from .shard import Shard, shard_region_string


@dataclass(frozen=True, slots=True)
class RunResult:
    success: bool
    output_path: Path
    elapsed_seconds: float
    returncode: int
    stderr: str
    timed_out: bool = False
    cached: bool = False


@dataclass(slots=True)
class RunnerConfig:
    ref_fasta: Path = field(default_factory=lambda: Path("testdata/hs37d5.fa"))
    bam_path: Path = field(default_factory=lambda: Path("testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"))
    java_jar: Path = field(default_factory=lambda: Path("VarDictJava/build/libs/VarDict-1.8.3.jar"))
    rust_bin: Path = field(default_factory=lambda: Path("target/debug-release/vardict"))
    base_dir: Path = field(default_factory=lambda: Path("tmp/na12878_parity"))
    parallel: int = 10
    java_parallel: int = 10
    java_heap: str = "8g"
    timeout: int = 600
    retry: int = 0
    freq: float = 0.01
    compress: bool = True
    rust_only: bool = False
    sample_name: str = "NA12878"


class ParityRunner:
    def __init__(self, project_root: Path, config: RunnerConfig):
        self.project_root = project_root
        self.config = config

    def run_java_shard(self, shard: Shard, test_config: TestConfig) -> RunResult:
        output_path = self._output_path("java", shard, test_config)
        if output_path.exists():
            return RunResult(True, output_path, 0.0, 0, "", cached=True)

        command = ["java"]
        if self.config.java_heap:
            command.append(f"-Xmx{self.config.java_heap}")
        command.extend(
            [
                "-classpath",
                str(self._resolve_path(self.config.java_jar)),
                "com.astrazeneca.vardict.Main",
                "-G",
                str(self._resolve_path(self.config.ref_fasta)),
                "-b",
                str(self._resolve_path(self.config.bam_path)),
                "-N",
                self.config.sample_name,
            ]
        )
        command.extend(self._freq_args(test_config))
        command.extend(["-th", "1"])
        command.extend(shlex.split(test_config.flags))
        command.extend(["-R", shard_region_string(shard)])
        return self._run_command(command, output_path=output_path, log_path=self._log_path("java", shard, test_config))

    def run_rust_shard(self, shard: Shard, test_config: TestConfig) -> RunResult:
        output_path = self._output_path("rust", shard, test_config)
        if output_path.exists() and not self.is_rust_stale(shard, test_config):
            return RunResult(True, output_path, 0.0, 0, "", cached=True)

        command = [
            str(self._resolve_path(self.config.rust_bin)),
            "-G",
            str(self._resolve_path(self.config.ref_fasta)),
            "-b",
            str(self._resolve_path(self.config.bam_path)),
            "-N",
            self.config.sample_name,
        ]
        command.extend(self._freq_args(test_config))
        command.extend(shlex.split(test_config.flags))
        command.extend(["-R", shard_region_string(shard)])
        return self._run_command(command, output_path=output_path, log_path=self._log_path("rust", shard, test_config))

    def is_java_cached(self, shard: Shard, test_config: TestConfig) -> bool:
        return self._output_path("java", shard, test_config).exists()

    def is_rust_stale(self, shard: Shard, test_config: TestConfig) -> bool:
        output_path = self._output_path("rust", shard, test_config)
        rust_bin = self._resolve_path(self.config.rust_bin)
        if not output_path.exists():
            return True
        if not rust_bin.exists():
            return True
        return rust_bin.stat().st_mtime > output_path.stat().st_mtime

    def _run_command(self, command: list[str], *, output_path: Path, log_path: Path) -> RunResult:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        attempts = self.config.retry + 1
        for attempt in range(1, attempts + 1):
            result = self._run_once(command, output_path=output_path, log_path=log_path)
            if result.success or attempt == attempts:
                return result
        return result

    def _run_once(self, command: list[str], *, output_path: Path, log_path: Path) -> RunResult:
        temp_handle = tempfile.NamedTemporaryFile(
            mode="wb",
            delete=False,
            dir=output_path.parent,
            prefix=f".{output_path.stem}.",
            suffix=".tmp",
        )
        temp_path = Path(temp_handle.name)
        temp_handle.close()

        start_time = time.monotonic()
        try:
            with temp_path.open("wb") as stdout_handle, log_path.open("wb") as stderr_handle:
                completed = subprocess.run(
                    command,
                    cwd=self.project_root,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    timeout=self.config.timeout,
                    check=False,
                )
            elapsed = time.monotonic() - start_time
        except subprocess.TimeoutExpired:
            elapsed = time.monotonic() - start_time
            temp_path.unlink(missing_ok=True)
            output_path.unlink(missing_ok=True)
            stderr_text = self._read_text(log_path)
            return RunResult(False, output_path, elapsed, 124, stderr_text, timed_out=True)

        stderr_text = self._read_text(log_path)
        if completed.returncode != 0:
            temp_path.unlink(missing_ok=True)
            output_path.unlink(missing_ok=True)
            return RunResult(False, output_path, elapsed, completed.returncode, stderr_text)

        if self.config.compress:
            with temp_path.open("rb") as source, gzip.open(output_path, "wb", compresslevel=1) as target:
                shutil.copyfileobj(source, target)
            temp_path.unlink(missing_ok=True)
        else:
            shutil.move(str(temp_path), str(output_path))
        return RunResult(True, output_path, elapsed, 0, stderr_text)

    def _freq_args(self, test_config: TestConfig) -> list[str]:
        tokens = shlex.split(test_config.flags)
        if "-f" in tokens:
            return []
        return ["-f", f"{self.config.freq}"]

    def _output_path(self, mode: str, shard: Shard, test_config: TestConfig) -> Path:
        suffix = ".tsv.gz" if self.config.compress else ".tsv"
        return self.config.base_dir / test_config.label / shard.chrom / mode / f"shard_{shard.label}{suffix}"

    def _log_path(self, mode: str, shard: Shard, test_config: TestConfig) -> Path:
        return self.config.base_dir / test_config.label / shard.chrom / mode / f"shard_{shard.label}.log"

    def _resolve_path(self, path: Path) -> Path:
        if path.is_absolute():
            return path
        return self.project_root / path

    @staticmethod
    def _read_text(path: Path) -> str:
        if not path.exists():
            return ""
        return path.read_text(encoding="utf-8", errors="replace")
