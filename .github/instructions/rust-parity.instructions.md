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
