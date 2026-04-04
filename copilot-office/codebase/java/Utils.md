# Utils

**Source**: `Utils.java`
**LOC**: 279
**Pipeline Stage**: Shared utility — called from every pipeline stage
**Rust counterpart**: `src/utils.rs` (partial), plus embedded copies in `src/mods/variant_realigner.rs` and `src/mods/structural_variants_processor.rs`
**Status**: complete

## Overview

`Utils` is a `final class` containing exclusively `static` utility methods. It provides:

1. **String join/format** — `toString()`, `join()`, `joinNotNull()`
2. **Perl-compatible string ops** — `substr()` (2 overloads), `charAt()` (2 overloads)
3. **Numeric formatting** — `roundHalfEven()`, `getRoundedValueToPrint()`
4. **Map utilities** — `getOrElse()` (with mutation side-effect)
5. **Parsing** — `toInt()`
6. **Sequence operations** — `complement()` (3 overloads), `reverse()`, `getReverseComplementedSequence()`
7. **Pattern matching** — `globalFind()`
8. **Collection** — `sum()`
9. **Error handling** — `printExceptionAndContinue()`

**Note**: `joinRef()`, `incCnt()`, `adjCnt()`, `findconseq()`, `strandBias()` are in `VariationUtils.java`, NOT in `Utils.java`.

## Method Inventory

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `toString(Collection)` | L19-L31 | yes | Space-delimited string from collection |
| `join(String, Object...)` | L39-L50 | yes | Delimiter-joined string from varargs |
| `joinNotNull(String, Object...)` | L58-L79 | yes | Join skipping nulls with special delimiter logic |
| `getOrElse(Map, K, V)` | L81-L87 | yes | Get-or-insert-default (**MUTATES map**) |
| `toInt(String)` | L89-L91 | yes | Wrapper for `Integer.parseInt()` |
| `roundHalfEven(String, double)` | L99-L103 | yes | Format→parse round-trip with HALF_EVEN |
| `getRoundedValueToPrint(String, double)` | L105-L109 | yes | Conditional formatting for output columns |
| `substr(String, int)` | L117-L123 | yes | Perl-compatible substring from index |
| `substr(String, int, int)` | L132-L148 | yes | Perl-compatible substring with length |
| `charAt(String, int)` | L155-L164 | yes | Perl-compatible char access, negative index |
| `charAt(StringBuilder, int)` | L172-L181 | yes | Same for StringBuilder |
| `sum(Collection)` | L188-L194 | yes | Sum of int-valued objects |
| `globalFind(Pattern, String)` | L201-L209 | yes | Collect all regex group(1) matches |
| `getReverseComplementedSequence(...)` | L217-L224 | yes | Reverse-complement a sub-read |
| `reverse(String)` | L226-L228 | yes | Reverse via StringBuffer |
| `complement(String)` | L230-L234 | yes | Per-base complement of string |
| `complement(byte[])` | L236-L241 | yes | In-place complement of byte array |
| `complement(char)` | L243-L247 | yes | Single character complement |
| `printExceptionAndContinue(...)` | L256-L279 | yes | Exception logging with count limit |

## Method Analyses

### 3.1 `toString(Collection<E>)` — L19-L31

**Signature**: `<E> String toString(Collection<E>)`

**Purpose**: Convert a collection to a space-delimited string.

**Algorithm**:
1. Get iterator. If empty, return `""`
2. Infinite loop: append element (or `"(this Collection)"` if element IS the collection to prevent recursion), break when done, append space between

**Edge cases**: Empty → `""`. Single element → no trailing space. Self-referencing collection → prints `"(this Collection)"`.

### 3.2 `join(String, Object...)` — L39-L50

**Signature**: `String join(String, Object...)`

**Purpose**: Join varargs with delimiter.

**Algorithm**: Loop args, append each, append delimiter between (not after last). Empty args → `""`.

**Parity note**: Java `StringBuilder.append(null_object)` produces literal `"null"`. If any caller passes null, output contains `"null"` string.

### 3.3 `joinNotNull(String, Object...)` — L58-L79

**Signature**: `String joinNotNull(String, Object...)`

**Purpose**: Join varargs, skipping nulls with nuanced delimiter logic.

**Algorithm**: For each arg:
- If `null`: append delimiter only if next is non-null; skip element
- If non-null: append element; append delimiter only if next is non-null

**Trace** of `joinNotNull("-", null, "A", null, "B")`:
- i=0: null, next="A" → append `"-"` → `"-"`
- i=1: "A", next=null → append "A", no delim → `"-A"`
- i=2: null, next="B" → append `"-"` → `"-A-"`
- i=3: "B", end → append "B" → `"-A-B"`

**Result**: `"-A-B"` — note the **leading delimiter** when first arg is null.

### 3.4 `getOrElse(Map<K,V>, K, V)` — L81-L87

**Signature**: `<K,V> V getOrElse(Map<K,V>, K, V)`

**Purpose**: Get value or insert default and return it.

```java
V v = map.get(key);
if (v == null) {
    v = or;
    map.put(key, v);  // ← MUTATES MAP
}
return v;
```

**CRITICAL**: Inserts default into map when key absent. Equivalent to Rust `entry().or_insert()`, NOT `get().unwrap_or()`.

**Null ambiguity**: `HashMap.get()` returns `null` for both "key absent" and "key present with null value". So if `map.put(key, null)` was called first, then `getOrElse` overwrites the null with the default.

### 3.5 `toInt(String)` — L89-L91

**Signature**: `int toInt(String)`

Direct `Integer.parseInt()` wrapper. Throws on null, empty, or non-numeric input.

### 3.6 `roundHalfEven(String, double)` — L99-L103

**Signature**: `double roundHalfEven(String, double)`

**Purpose**: Round a double using HALF_EVEN (banker's rounding) via format→parse round-trip.

```java
return Double.parseDouble(new DecimalFormat(pattern).format(value));
```

**Step-by-step**:
1. `DecimalFormat(pattern)` — default `RoundingMode.HALF_EVEN`
2. `format(value)` → string at pattern's precision
3. `Double.parseDouble(result)` → re-encode as double

**Patterns used in codebase**:

| Pattern | Decimals | Used For |
|---------|----------|----------|
| `"0.0"` | 1 | Quality, mapq, position mean, NM |
| `"0.000"` | 3 | MSI, q-ratio |
| `"0.0000"` | 4 | Frequency, extra freq, HQ reads freq |
| `"0.00000"` | 5 | P-values |

**HALF_EVEN examples**:
- `roundHalfEven("0.0", 2.25)` → `"2.2"` → `2.2` (tie → even digit 2)
- `roundHalfEven("0.0", 2.35)` → `"2.4"` → `2.4` (tie → even digit 4)
- `roundHalfEven("0.0", 2.45)` → `"2.4"` → `2.4` (tie → even digit 4)

The format→parse round-trip means the result carries double-precision representation artifacts from the parse step.

### 3.7 `getRoundedValueToPrint(String, double)` — L105-L109

**Signature**: `String getRoundedValueToPrint(String, double)`

**Purpose**: Format a double for TSV output with conditional formatting.

```java
return value == Math.round(value)
    ? new DecimalFormat("0").format(value)          // whole number → "42"
    : new DecimalFormat(pattern).format(value)
        .replaceAll("0+$", "");                     // fractional → strip trailing zeros
```

**Branches**:
- `value == Math.round(value)` ← whole number test. `Math.round(double)` returns `long`.
  - True: format with `"0"` pattern → `"0"`, `"1"`, `"42"` (no decimal point)
  - False: format with pattern, then strip trailing `0` chars via regex

**Edge cases**:
- `value = 0.0` → `0.0 == 0L` → true → `"0"`
- `value = 1.0` → true → `"1"`
- `value = 1.5` → `Math.round(1.5) == 2L`, `1.5 != 2` → false
- `value = 0.1230` with `"0.0000"` → `"0.1230"` → strip → `"0.123"`
- `value = 0.1000` with `"0.0000"` → `"0.1000"` → strip → `"0.1"`
- **Trailing dot edge**: if `value ≈ 1.0000001`, `Math.round()` gives 1, `1.0000001 != 1` → false; `DecimalFormat("0.0").format()` → `"1.0"`; `replaceAll("0+$","")` → `"1."` — **trailing decimal point**

### 3.8 `substr(String, int)` — L117-L123

**Signature**: `String substr(String, int)`

**Purpose**: Perl-compatible `substr($str, $idx)`.

**Algorithm**:
- `idx >= 0`: `string.substring(min(length, idx))` — clamped, never throws
- `idx < 0`: `string.substring(max(0, length + idx))` — counts from right, clamped to 0

**Examples**:

| Call | Result | Why |
|------|--------|-----|
| `substr("ABCDE", 2)` | `"CDE"` | From index 2 |
| `substr("ABCDE", 10)` | `""` | Clamped to length |
| `substr("ABCDE", -2)` | `"DE"` | From position 3 |
| `substr("ABCDE", -10)` | `"ABCDE"` | Clamped to 0 |
| `substr("", 0)` | `""` | |
| `substr("", -1)` | `""` | max(0, -1) → 0 |

### 3.9 `substr(String, int, int)` — L132-L148

**Signature**: `String substr(String, int, int)`

**Purpose**: Perl-compatible `substr($str, $begin, $len)`.

**Algorithm**:
1. If `begin < 0`: `begin = length + begin` (**NO CLAMPING** — can go negative → exception)
2. If `len > 0`: `substring(begin, min(begin+len, length))` — end clamped
3. If `len == 0`: return `""`
4. If `len < 0`: `end = length + len`; if `end < begin` return `""`; else `substring(begin, end)`

**Critical difference from 1-arg**: The 1-arg overload clamps negative begin with `Math.max(0, ...)`. The 2-arg overload does NOT — negative begin where `|begin| > length` produces a negative `begin` and `StringIndexOutOfBoundsException`.

**Examples**:

| Call | Result | Why |
|------|--------|-----|
| `substr("ABCDE", 1, 3)` | `"BCD"` | begin=1, end=4 |
| `substr("ABCDE", 0, 10)` | `"ABCDE"` | end clamped to 5 |
| `substr("ABCDE", -2, 1)` | `"D"` | begin=3, len=1 |
| `substr("ABCDE", -2, 2)` | `"DE"` | begin=3, len=2 |
| `substr("ABCDE", 1, -1)` | `"BCD"` | begin=1, end=4 |
| `substr("ABCDE", 3, -1)` | `"D"` | begin=3, end=4 |
| `substr("ABCDE", 1, 0)` | `""` | len==0 |
| `substr("ABCDE", -10, 2)` | **THROWS** | begin=-5, neg index |

### 3.10 `charAt(String, int)` — L155-L164

**Signature**: `char charAt(String, int)`

**Purpose**: Perl-compatible char access with negative indexing.

**Algorithm**:
- `index < 0`: `i = length + index`. If `i < 0` → return `(char)-1` (`\uFFFF`). Else `str.charAt(i)`.
- `index >= 0`: `str.charAt(index)` — throws if out-of-bounds.

**Critical**: Negative overflow → sentinel `(char)-1`. Positive overflow → exception. This asymmetry is intentional.

| Call | Result |
|------|--------|
| `charAt("ABCDE", 0)` | `'A'` |
| `charAt("ABCDE", -1)` | `'E'` |
| `charAt("ABCDE", -5)` | `'A'` |
| `charAt("ABCDE", -6)` | `(char)-1` = `\uFFFF` |
| `charAt("ABCDE", 5)` | **THROWS** |

### 3.11 `charAt(StringBuilder, int)` — L172-L181

**Signature**: `char charAt(StringBuilder, int)`

Identical logic to 3.10 but for `StringBuilder`.

### 3.12 `sum(Collection<?>)` — L188-L194

**Signature**: `int sum(Collection<?>)`

Sum of int-valued objects: `result += toInt(String.valueOf(object))`. Empty → 0.

### 3.13 `globalFind(Pattern, String)` — L201-L209

**Signature**: `List<String> globalFind(Pattern, String)`

Collects all `group(1)` matches from a jregex Pattern. Returns `LinkedList<String>`. Note: uses **jregex**, not java.util.regex.

### 3.14 `getReverseComplementedSequence(SAMRecord, int, int)` — L217-L224

**Signature**: `String getReverseComplementedSequence(SAMRecord, int, int)`

Extracts `[startIndex, startIndex+length)` from read bases, reverse-complements in place via `SequenceUtil.reverseComplement()`. Negative `startIndex` counts from end.

### 3.15 `reverse(String)` — L226-L228

**Signature**: `String reverse(String)`

`new StringBuffer(string).reverse().toString()` — simple byte reversal for ASCII genomic data.

### 3.16-3.18 `complement()` overloads — L230-L247

All three overloads delegate to `SequenceUtil.complement(byte)`:

| Input | Output |
|-------|--------|
| `A` | `T` |
| `T` | `A` |
| `C` | `G` |
| `G` | `C` |
| `a` | `t` |
| `t` | `a` |
| `c` | `g` |
| `g` | `c` |
| **anything else** | **unchanged** (pass-through) |

This includes `N` → `N`, digits → digits, etc.

### 3.19 `printExceptionAndContinue(...)` — L256-L279

**Signature**: `void printExceptionAndContinue(Exception, String, String, Region)`

Logs exception to stderr, increments atomic counter, throws if `> MAX_EXCEPTION_COUNT`. Two message formats depending on whether `region` is null.

## Known Parity Traps

### Trap 1: `substr()` 2-arg negative begin has NO clamping

`substr(str, begin, len)` with `begin < 0` computes `begin = length + begin` without `Math.max(0, ...)`. If `|begin| > length`, this goes negative and throws. The 1-arg overload DOES clamp. This asymmetry must be preserved.

### Trap 2: `charAt()` returns `(char)-1` for negative overflow

`charAt(str, -N)` where `N > length` returns `\uFFFF` (65535). In Rust with `u8` bytes, there's no equivalent. The Rust port uses `Option<u8>` → `None`. Callers must handle `None` the same way Java handles `\uFFFF`.

### Trap 3: `charAt()` positive OOB throws, negative OOB returns sentinel

Only negative underflow returns `(char)-1`. Positive overflows throw `StringIndexOutOfBoundsException`. The Rust port must distinguish these two failure modes.

### Trap 4: `roundHalfEven()` format→parse double-conversion

The format→parse round-trip means results carry double-precision artifacts. `0.00005` in double is slightly above 0.00005, so HALF_EVEN rounds up instead of to-even. The Rust implementation must replicate this exactly.

### Trap 5: `roundHalfEven()` pattern determines precision

Pattern `"0.0"` = 1 decimal, `"0.0000"` = 4 decimals. The Rust implementation parses the pattern string to extract scale.

### Trap 6: `getRoundedValueToPrint()` trailing-zero stripping

`replaceAll("0+$", "")` strips all trailing zeros. Can produce trailing decimal point (e.g., `"1."`) when a value very close to a whole number rounds to exactly that whole number but fails the `value == Math.round(value)` test.

### Trap 7: `getOrElse()` mutates the map

Inserts default when key absent. Changes map size and `LinkedHashMap` iteration order. Use Rust `entry().or_insert()`.

### Trap 8: `getOrElse()` null-value ambiguity

`HashMap.get()` returns `null` for both "key absent" and "key has null value". `getOrElse` cannot distinguish them — it overwrites null-valued entries with the default.

### Trap 9: `complement()` preserves case and passes through non-ACGT unchanged

htsjdk `SequenceUtil.complement()`: `A↔T`, `C↔G` (preserving case). All other bytes pass through **unchanged** — including `N`, digits, special chars.

**Current Rust parity break**: `complement_base_u8` in `structural_variants_processor.rs` does `to_ascii_uppercase()` (loses case) and maps unknown → `b'N'` (should pass through unchanged).

### Trap 10: `joinNotNull()` leading delimiter on null-first

When first arg is null and second is non-null, a delimiter appears at the start of output.

### Trap 11: `substr()` 2-arg negative len end calculation

`len < 0`: `end = length + len`. Not clamped. Combined with negative-begin conversion, both transformations interact. Must trace at each call site.

### Trap 12: `globalFind()` uses jregex, not java.util.regex

Different regex engine. Simple patterns work the same; complex features may differ.

### Trap 13: `join()` with null elements produces literal `"null"`

Java `StringBuilder.append(null_object)` → `"null"`.

### Trap 14: `reverse()` uses StringBuffer

Functionally equivalent to StringBuilder for ASCII data. No parity concern but documents thread safety intent.

## Cross-Module Dependencies

**Called by (extensively)**:
- **CigarParser**: `substr()`, `charAt()`, `complement()`, `reverse()`, `roundHalfEven()`
- **VariationRealigner**: `substr()`, `charAt()`, `complement()`, `reverse()`, `join()`, `globalFind()`
- **StructuralVariantsProcessor**: `substr()`, `charAt()`, `complement()`, `reverse()`, `getReverseComplementedSequence()`
- **ToVarsBuilder**: `roundHalfEven()`, `getOrElse()`, `substr()`, `charAt()`
- **OutputVariant printers**: `getRoundedValueToPrint()`, `join()`, `joinNotNull()`
- **SAMFileParser/Modes**: `printExceptionAndContinue()`

**Calls externally**:
- `htsjdk.samtools.util.SequenceUtil.complement(byte)` — base complement
- `htsjdk.samtools.util.SequenceUtil.reverseComplement(byte[])` — reverse + complement
- `htsjdk.samtools.util.StringUtil` — string↔byte conversions
- `jregex.Pattern` / `jregex.Matcher` — regex engine
- `java.text.DecimalFormat` — number formatting (HALF_EVEN default)
- `Configuration.MAX_EXCEPTION_COUNT` — exception limit
- `GlobalReadOnlyScope.instance().conf` — global configuration

## Float Formatting Summary

| Pattern | Decimals | Used For | Rounding |
|---------|----------|----------|----------|
| `"0"` | 0 | Whole numbers in `getRoundedValueToPrint` | N/A |
| `"0.0"` | 1 | Quality, mapq, position mean, NM | HALF_EVEN |
| `"0.000"` | 3 | MSI, q-ratio | HALF_EVEN |
| `"0.0000"` | 4 | Frequency, extra freq, HQ reads freq | HALF_EVEN |
| `"0.00000"` | 5 | P-values | HALF_EVEN |
