---
description: "Use when updating tests, fixture snapshots, parity cases, JSONL dumps, or output-affecting Rust code that needs validation. Covers narrowest-test selection, parity manifest usage, and debug artifact workflow."
name: "Parity And Fixture Testing"
applyTo: "tests/**"
---

# Parity And Fixture Testing

- Start with the narrowest useful validation for the code path you changed, then widen only if needed.
- Prefer stage-level fixture tests in `tests/*_fixture_test.rs` before full parity runs in `tests/integration_test.rs`.
- Use `tests/parity_case_manifest.csv` to select `RUN_NOW` cases and relevant tags before running integration parity checks.
- When a change affects output ordering, formatting, numeric rendering, CLI compatibility, or cross-stage pipeline behavior, include a parity-oriented validation step rather than only unit coverage.
- Use the existing JSONL snapshot workflow for pipeline debugging instead of ad hoc print debugging.
- Keep debug artifacts under `./tmp/` and use the existing dump environment variables when you need intermediate stage output: `VARDICT_CIGAR_PARSER_JSONL`, `VARDICT_TO_VARS_JSONL`, `VARDICT_VARIANT_REALIGNER_JSONL`, and `VARDICT_STRUCTURAL_VARIANTS_JSONL`.
- Useful parity test controls: `VARDICT_RUN_NOW_SIMPLE_LIMIT`, `VARDICT_RUN_NOW_FULL_SWEEP`, `VARDICT_RUN_NOW_STRICT`, `VARDICT_TEST_LOG`, and `TIER1_DUMP_LINES`.
- Module-to-fixture mapping is direct: `cigar_parser` -> `cigar_parser_fixture_test`, `to_vars_builder` -> `tovars_fixture_test`, `variant_realigner` -> `variant_realigner_fixture_test`, `structural_variants_processor` -> `structural_variants_fixture_test`.
- Prefer `cargo test --profile debug-release` for validation work unless you are explicitly checking optimized behavior.