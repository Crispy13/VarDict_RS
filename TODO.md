
## [ ] Fix: Output Storing Behavior of Rust (written on 2026-03-17)
Good question. Here's my assessment:

**The batching fix is a workaround, not a proper fix.** It bounds memory at the *dispatch* layer by chunking regions, but the core pipeline in `process_regions_vardict()` still accumulates all results within each batch into `Vec<RegionResult>` before printing. Java doesn't work this way — it uses a bounded `BlockingQueue(10)` that streams results to stdout as they complete, holding at most ~10 results in memory at any time.

**For the test script scenario, the current fix is sufficient.** Each shard invocation uses `-R region` (a single region), so there's only 1 region per batch. No accumulation problem.

**For real-world BED file usage, the batching helps but isn't ideal.** With 4 threads, batch size = 16. That's 16 `RegionResult`s (each holding `Vec<String>` of output lines) accumulated before any are printed. Manageable, but still fundamentally different from Java's streaming model.

**The proper fix** would be to replace the collect-then-print pattern with a bounded channel in `parallel_pipeline.rs`:
- Producer threads (rayon workers) push `RegionResult` into a `crossbeam::channel::bounded(N)` 
- A dedicated consumer thread prints and drops results as they arrive
- Memory stays at O(channel_capacity), matching Java's `BlockingQueue` semantics
- The batching workaround in `select_region_batches_for_execution` becomes unnecessary

**My recommendation**: The current fix unblocks the immediate problem (test suite memory explosion). The streaming fix is worth doing, but it's a separate concern — it affects the production binary's behavior with large BED files, not the test infrastructure. I'd prioritize it after the Tier 1 parity bugs are resolved, since those are the gate to correctness validation. The streaming fix is a performance/robustness improvement that doesn't affect output correctness.

Should I add it to the project plan as a follow-up task?

Answered the question about Rust's result storing behavior: the current batching fix is a workaround that's sufficient for the test suite, but the proper fix would be a bounded channel in `parallel_pipeline.rs` to match Java's streaming `BlockingQueue` model. Recommended deferring the streaming fix until after Tier 1 parity bugs are resolved. 

