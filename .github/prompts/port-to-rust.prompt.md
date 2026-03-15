---
description: "Port a VarDictJava method to Rust with output parity. Translates Java algorithm logic to idiomatic Rust while maintaining byte-identical output behavior."
# agent: "rust-implementer"
argument-hint: "Specify the Java method to port, e.g. CigarParser.parseCigar, or describe the parity fix needed"
---

Port the specified VarDictJava method to Rust, or fix an existing Rust implementation for parity.


1. **Read the existing Rust code** in the target module first to understand current structure
2. **Follow the type mapping** from `rust-parity.instructions.md`:
   - `int` → `i32`, `long` → `i64`, `double` → `f64`
   - `LinkedHashMap` → `IndexMap`, `TreeMap` → `BTreeMap`
   - `null` → `Option::None`
3. **Preserve every branch** from the Java control flow — do not optimize away any path
4. **Add traceability**: Document which Java class and method this ports, with line numbers
5. **Match float formatting** exactly (Java `DecimalFormat` behavior)
6. **Use `wrapping_*` arithmetic** where Java integers could overflow
7. **Run `cargo check`** after implementation to verify compilation

Follow `rust.instructions.md` for general Rust style and `rust-parity.instructions.md` for parity-specific rules.

When parity correctness conflicts with idiomatic Rust, choose parity and add a comment explaining why.
