[ ] MAP-001: Align remaining CLI argument surface with Java `CmdParser`
    - Java authority: `VarDictJava/src/main/java/com/astrazeneca/vardict/CmdParser.java`
    - Rust targets: `src/bin/vardict.rs`, `src/conf.rs`
    - Validation: `cargo test --bin vardict`

[ ] MAP-002: Replace residual hard-coded constants with `Configuration` fields
    - Java authority: `VarDictJava/src/main/java/com/astrazeneca/vardict/Configuration.java`
    - Rust target: `src/conf.rs`
    - Validation: `cargo test --lib conf`

[ ] MAP-003: Keep integration option parsing parity for compact/legacy aliases
    - Java authority: `VarDictJava/src/test/java/com/astrazeneca/vardict/integrationtests/IntegrationTest.java`
    - Rust target: `tests/integration_test.rs`
    - Validation: `cargo test --test integration_test test_manifest_tier1_simple_raw_rust_vs_java_first_mismatch -- --ignored --nocapture`

Reference queue doc: `copilot-desk/java-rust-mapping-execution-2026-02-24.md`
