# Benchmarks

This document summarizes a focused end-to-end benchmark of VarDict-rs against VarDictJava on the same dataset, reference, genomic region, and command-line options.

## Methodology

- Workload: a 1 MB region on chromosome 1, 1:100000000-101000000.
- Dataset: NA12878 low-coverage whole-genome sequencing BAM against the hs37d5 reference.
- Modes: standard variant calling and pileup.
- Iterations: 1 warmup run per implementation and mode, followed by 5 timed iterations.
- Concurrency: single-threaded for both implementations.
- Timing: wall-clock runtime measured with the bash time command, using the real field.
- Build and runtime settings: Rust used the debug-release profile; Java used OpenJDK 8 with an 8 GB heap limit.
- Output handling: benchmark runs wrote output to /dev/null so the measurements reflect execution time rather than terminal or file I/O.

## Environment

| Item | Value |
| --- | --- |
| CPU | AMD Ryzen 7 5800X 8-Core Processor |
| Visible cores / threads | 6 / 12 in WSL2 |
| Memory | 24 GB |
| OS | Ubuntu 24.04.4 LTS on WSL2, kernel 6.6.87.2-microsoft-standard-WSL2 |
| Rust | 1.94.0 (2026-03-02) |
| Java | OpenJDK 1.8.0_482 |
| Build profile | debug-release (release optimizations with debug symbols) |
| Dataset | NA12878 low-coverage WGS, 16 GB BAM |
| Reference | hs37d5 |

## Results

| Workload | Rust mean ± stddev (s) | Rust min (s) | Rust max (s) | Java mean ± stddev (s) | Java min (s) | Java max (s) | Speedup (Java / Rust) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Variant calling | 3.17 ± 0.22 | 2.98 | 3.60 | 5.60 ± 0.03 | 5.56 | 5.66 | 1.77x |
| Pileup | 7.94 ± 0.48 | 7.29 | 8.60 | 15.37 ± 0.58 | 14.44 | 16.27 | 1.94x |

## Output Parity

For the benchmark commands above, VarDict-rs and VarDictJava produced byte-identical output for the same inputs. Broader parity validation across additional datasets and option combinations is tracked separately in the repository parity workflow.

Results are hardware-dependent and should be treated as end-to-end measurements for this specific configuration.