---
description: Cross-cutting operational policies for VarDict-rs work.
applyTo: '**'
---

# Operational Policies

## Build Profile
- Use `cargo build --profile debug-release` for development, debugging, and verification.
- Reserve `cargo build --release` for production or explicit release validation.

## Environment
- Activate `rust_build_env` before builds or tests, then set `LIBCLANG_PATH=$CONDA_PREFIX/lib`.
- If the conda env is broken, stop and ask whether to continue with Conda or switch to a Python `venv`.

## Temporary Files
- Create all temporary or intermediate files under `./tmp`.
- Do not use the system `/tmp` directory.

## Test Command
- Use `cargo test --profile debug-release -- --include-ignored` for validation.
- Include ignored tests because parity regressions are often parked there first.

## Lint
- Run `cargo clippy -- -D warnings` and `cargo fmt --check` before closing a Rust change.

## Sweep Policy
- Never use `--no-stop` for multi-config sweeps; stop on the first real parity mismatch, fix it, then re-run.
- `--no-stop` is acceptable only for single-config, single-chromosome runs or when the failure is confirmed to be non-Rust, for example Java OOM.

## Tiered Config Order
- Run cheaper, high-signal configs first: core pipeline, then filter or variant options, then combo configs, then no-realign or debug, then pileup last.

## Commit Discipline
- Land the fix, Java-derived fixture, and regression test in the same commit.

## Skill Discovery
- Repository skills live under `.github/skills/{name}/SKILL.md`.
