---
description: Preferences for Rust projects in addition to awesome copilot directions.
applyTo: '**/*.rs'
---
1. Do not modify Cargo.lock file directly. It should be modified by cargo commands.

2. Prefer `tracing` and `tracing-subscriber` for logging over `println!` or `eprintln!` for better performance and flexibility.  
  
Make at least two subscriber init functions:
```rust
# one for development including callsite. (e.g. test)
.with_timer(ChronoLocal::rfc_3339())
.with_file(true)
.with_line_number(true)
.with_target(true)
.with_ansi(false)

# one for production without callsite. (CLI normal behavior)
.with_timer(ChronoLocal::rfc_3339())
.with_file(false)
.with_line_number(false)
.with_target(false)
```

3. Use `debug-release` profile for testing and debugging. If Cargo.toml does not have it, ask user to add it.

4. Prefer enums or generic code over dynamic dispatch (trait objects) for better performance and compile-time checks.