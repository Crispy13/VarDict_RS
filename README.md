# VarDictRust

VarDictRust is a Rust port of VarDictJava targeting improved performance while giving the same output.

This repository is being released as a public beta. It is useful for experimentation, benchmarking, and parity work, but it is not yet 100% byte-identical to VarDictJava.

This repository has been written heavily by VS Code Copilot, with human review and editing.

## Quick start

1. Install Rust
2. Build: `cargo build --profile deploy`
3. Run: `./target/deploy/vardict`

## Note

- Output parity is still in progress.
- Some option combinations and edge cases can still diverge from Java.
- Treat the project as beta/experimental for now.
- For critical comparisons, validate Rust output against VarDictJava.

## Current parity progress

- Testing the NA12878 mapped low-coverage BAM against the hs37d5 reference as the main parity dataset.
- Parity sweeps are still in progress, so the repository should be treated as a public beta.
- See [Parity workflow](docs/PARITY_WORKFLOW.md) for the current test flow and BAM status.
