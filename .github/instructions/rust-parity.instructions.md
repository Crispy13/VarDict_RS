---
description: "Use when writing or modifying Rust code for VarDictJava parity, implementing ported methods, fixing parity mismatches, or translating Java algorithms to Rust. Covers type mapping, float formatting, collection ordering, and output-matching techniques."
applyTo: "**/*.rs"
---

# VarDictJava Parity — Rust Implementation Rules

These rules supplement `rust.instructions.md` with parity-specific requirements. When in conflict, **parity correctness takes precedence** over idiomatic Rust style.

## Type Mapping (Java → Rust)

| Java Type | Rust Type | Notes |
|-----------|-----------|-------|
| `int` | `i32` | Match signed 32-bit semantics |
| `long` | `i64` | |
| `double` | `f64` | |
| `float` | `f32` | Rare in VarDict |
| `boolean` | `bool` | |
| `char` | `u8` | Genomic data is ASCII; use `u8` not `char` for base operations |
| `String` | `String` / `&str` | |
| `byte[]` | `Vec<u8>` / `&[u8]` | |
| `int[]` | `Vec<i32>` / `&[i32]` | |
| `HashMap<K,V>` | `HashMap<K,V>` | Only when order doesn't matter |
| `LinkedHashMap<K,V>` | `IndexMap<K,V>` | **Critical**: preserves insertion order |
| `TreeMap<K,V>` | `BTreeMap<K,V>` | |
| `ArrayList<T>` | `Vec<T>` | |
| `null` | `Option::None` | Every Java null check → `Option` match |
| `StringBuilder` | `String` with `push_str` | |

## Float Formatting (Critical for Parity)

Java's `DecimalFormat("0.0000")` uses **HALF_EVEN** (banker's) rounding. Rust's default `format!("{:.4}", v)` uses **HALF_TO_EVEN** which matches. Verify with edge cases:

```rust
// Java: new DecimalFormat("0.0000").format(0.00005) → "0.0000"  (rounds to even)
// Rust: format!("{:.4}", 0.00005_f64) → "0.0001"  ← MISMATCH!
// Solution: implement Java-compatible rounding when needed
```

When exact match fails, implement a `java_format` helper that replicates `DecimalFormat` behavior:
```rust
fn java_format_double(value: f64, decimal_places: usize) -> String;
```

## Integer Overflow

Java silently wraps on overflow. Rust panics in debug mode. For arithmetic that could overflow in Java:
```rust
// Use wrapping arithmetic to match Java behavior
let result = a.wrapping_add(b);
let result = a.wrapping_mul(b);
```

## Collection Ordering Rules

- **Output-affecting maps**: MUST use `IndexMap` (insertion-ordered)
- **Internal-only maps**: `HashMap` is acceptable if iteration order never reaches output
- **Sorted maps**: Use `BTreeMap` for Java `TreeMap`
- **When unsure**: Use `IndexMap` — correctness over performance

## Traceability

Every Rust function that ports a Java method must include:
```rust
/// Ported from: `com.astrazeneca.vardict.modules.CigarParser.parseCigar()`
/// Java source: CigarParser.java:L142-L380
fn parse_cigar(...) { ... }
```

## Parity-Critical Patterns

### Variant Description String Handling
Maintain the exact same string encoding as Java:
- `+SEQ` for insertions, `-N` for deletions, `#SEQ` for complex
- String concatenation order must match
- Regex patterns must produce identical matches

### Position Arithmetic
Genomic position calculations are 0-based or 1-based depending on context:
- BED coordinates: 0-based half-open
- SAM/BAM positions: 1-based
- VarDict internal: mixed — follow Java exactly

### Conditional Logic Preservation
When Java has nested `if/else` chains, preserve the exact branching structure:
```rust
// DO: Match Java's branch order
if condition_a {
    // Java branch A logic
} else if condition_b {
    // Java branch B logic
} else {
    // Java default logic
}

// DON'T: Refactor into match/early-return unless proven equivalent
```

### Output Column Order
Tab-delimited output columns must match exactly:
- Simple mode: 36 columns
- Amplicon mode: 38 columns
- Somatic mode: 55 columns

Verify column count and order against Java `*OutputVariant` classes.

## Debugging Parity Mismatches

When output differs from Java:
1. Identify the first differing column in the output line
2. Trace that column back to the Java `OutputVariant` printer
3. Find which `Variant` field populates it
4. Trace the field back through `ToVarsBuilder` → `CigarParser`
5. Compare the Rust equivalent at each stage
6. Check: float formatting? collection ordering? off-by-one? null handling?

## Find → Fix → Test Rule

**One parity failure = one new named `#[test]` function.** Every parity bug fix must include a regression test that locks in the fix. The workflow:

1. **Find**: Identify the output diff (expected Java vs actual Rust)
2. **Minimize**: Reduce to the smallest input (region, reads, options) that reproduces the diff
3. **Fixture**: Extract the expected output **from Java** and save it as the test fixture (`.tsv` / `.tsv.gz` under `tests/fixtures/` or inline in the test). **Never generate reference output from Rust** — that locks in Rust bugs instead of catching them.
4. **Test first**: Add exactly one new `#[test]` function in `tests/integration_test.rs` (or the relevant fixture test file) that compares Rust output against the Java fixture. The test **must fail before the fix** and pass after it. Requirements:
   - Is `#[ignore]`d with a descriptive reason string
   - Is named following this convention:
     - For NA12878 BAM parity: `test_target_bam_{bam_slug}_{chr}_{description}_parity`
       - Example: `test_target_bam_na12878_low_coverage_chr11_raw_parity`
     - For integration testcase files: `test_parity_{mode}_{case_slug}`
     - For unit-level module bugs: `test_{module}_{description}_parity`
       - Example: `test_cigar_parser_insertion_at_position_42_parity`
5. **Fix**: Correct the Rust code to match Java behavior — the test from step 4 is your acceptance gate
6. **Verify**: Run `cargo test --profile debug-release -- --include-ignored` to confirm the new test passes and no existing tests broke
7. **Commit**: Fix, fixture, and test go in the same commit

This prevents regressions where fixing one parity issue silently re-breaks another.

See `.github/instructions/ops-policy.instructions.md` for sweep strategy and tiered config test order.

## Codebase Cache

A progressive Rust codebase cache is maintained at `copilot-office/codebase/rust/VarDict-rs-CODEBASE.md`. Before deeply analyzing a Rust module:

1. **Check the cache first**: Read the relevant module file linked from the index. If architecture, Java correspondence, and parity traps are documented, use them to orient quickly. Verify key facts against actual source if critical.
2. **Analyze from source if missing**: If the module file doesn't exist or lacks the information you need, analyze the Rust source directly.
3. **Update the cache after significant work**: If you gained architectural insight (new parity traps, divergences, dependency maps), update the module's cache file. This is a secondary task — always complete your primary work first.
