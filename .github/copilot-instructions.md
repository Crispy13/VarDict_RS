# Project Guidelines

## Mission

- This repository ports VarDictJava to Rust.
- The default objective for code changes is exact Java output parity first, then better runtime than Java.
- Do not trade output parity for speed unless the change is explicitly marked, measured, and validated against parity tests.

## Architecture

- Keep the main pipeline model in mind: BAM/read processing flows through `cigar_parser` -> `to_vars_builder` -> realignment and structural-variant processing -> output formatting.
- `src/bin/vardict.rs` owns CLI compatibility and argument handling. Reusable behavior belongs in library modules under `src/mods/`, `src/data/`, `src/variants/`, and `src/scopedata/`.
- `src/mods/vardict_pipeline.rs` and `src/mods/pipeline.rs` coordinate end-to-end behavior. Changes there can affect ordering, formatting, and parity across many tests.
- `src/scopedata/global_read_only_scope.rs` and `src/data/shared_reference.rs` provide process-wide state and cached reference data. Treat them as shared infrastructure, not local implementation detail.

## Build And Test

- Use the `rust_build_env` conda environment for repo work. If conda activation fails under strict shell settings, surface that clearly instead of guessing.
- Prefer `cargo test --profile debug-release` for validation and debugging. Use `cargo build --release` or `cargo build --profile deploy` only when checking optimized behavior.
- Start with the narrowest useful test scope: fixture tests in `tests/*_fixture_test.rs` for stage-level behavior, then parity coverage in `tests/integration_test.rs`.
- Use `tests/parity_case_manifest.csv` to understand which Java-backed parity cases are expected to run now and which tags they cover.
- Keep temporary outputs and debug artifacts under `./tmp/`, never system `/tmp`.

## Conventions

- Follow the more specific rules in `.github/instructions/*.md`; keep this file limited to repo-wide guidance.
- Preserve Java-compatible output details, including record ordering, formatting, and numeric text rendering. Exact text parity matters in this project.
- Use `tracing`-based logging, not `println!`, for new diagnostic output.
- Prefer measuring before optimizing. Check repository memory notes before retrying known performance ideas that already regressed or proved unstable.
- When debugging pipeline stages, prefer the existing JSONL fixture and snapshot workflow over ad hoc print-based debugging.