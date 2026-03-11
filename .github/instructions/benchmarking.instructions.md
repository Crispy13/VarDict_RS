---
description: "Use when optimizing hot paths, benchmarking runtime, comparing Rust against Java VarDict, or changing hashers, allocators, parallelism, or batching. Covers benchmark hygiene, representative workloads, and known repo-specific pitfalls."
name: "Benchmarking And Optimization"
---

# Benchmarking And Optimization

- Preserve exact Java output parity before pursuing speed. If parity is uncertain, validate first and state the risk explicitly.
- Measure a before-and-after baseline and use the same workload, profile, and environment for both runs.
- Prefer representative workloads over tiny synthetic shards when making decisions about performance changes.
- Treat threading, batching, and allocator changes as empirical. Do not assume they help in this codebase without measurement.
- Keep benchmark outputs, temporary summaries, and ad hoc scripts under `./tmp/`.
- Confirm long-running benchmark jobs actually finished before trusting summary files or reported timings.
- Consult repository memory before retrying known regressions or unstable ideas, especially around global hasher changes, variant realigner deletion lookup shortcuts, and `mimalloc-global`.
- Use `cargo test --profile debug-release` for validation and debugging, then move to optimized builds only when the behavior is already understood.
- Call out when a performance change can affect record ordering, formatting, or any parity-sensitive behavior, even if the runtime result looks better.