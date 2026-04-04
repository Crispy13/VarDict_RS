# output_variant

**Source**: `src/mods/output_variant.rs`
**LOC**: ~2,677
**Java counterpart**: `SimpleOutputVariant.java`, `SomaticOutputVariant.java`, `AmpliconOutputVariant.java` → [Java cache](../java/OutputVariant.md)
**Status**: complete
**Last verified**: 2026-04-04

## Overview

Final output formatting stage of the VarDict pipeline. Converts in-memory `Variant` objects into tab-delimited text lines for three modes: **Simple** (36/38 columns), **Somatic** (55/61 columns), and **Amplicon** (38/40 columns). Column counts depend on whether fisher-mode is enabled. Also includes a full `FisherExact` test implementation (ported from Apache Commons Math via Java's `FisherExact.java`).

## Public API

| Symbol | Purpose |
|--------|---------|
| `SimpleOutputVariant` | 36/38-column single-sample output struct |
| `AmpliconOutputVariant` | 38/40-column amplicon output struct |
| `SomaticOutputVariant` | 55/61-column tumor/normal output struct |
| `SimpleOutputVariant::from_variant()` | Constructor from `Variant` |
| `SimpleOutputVariant::empty()` | Placeholder sentinel row (genotype="0", bias="0") |
| `SimpleOutputVariant::empty_null_variant()` | Null-variant row (genotype="", bias="0;0") |
| `to_string_36_columns()` / `to_string_38_columns()` | Non-fisher / fisher formatters |
| `AmpliconOutputVariant::from_variant()` | Constructor with amplicon metadata |
| `SomaticOutputVariant::from_variants()` | Constructor from 4 optional Variants |
| `get_*_header_line()` | Static header strings per mode |
| `FisherExact` | Hypergeometric p-value and odds ratio |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `SimpleOutputVariant` | `SimpleOutputVariant.java` | 1:1 port |
| `AmpliconOutputVariant` | `AmpliconOutputVariant.java` | 1:1; region logic differs slightly |
| `SomaticOutputVariant` | `SomaticOutputVariant.java` | 1:1; constructor takes 4 Variants |
| `FisherExact` | `FisherExact.java` + Apache Commons Math | Full port of hypergeometric CDF |
| `format_rounded_value_to_print()` | `Utils.getRoundedValueToPrint()` | Fisher-mode formatting |
| `format_f64()` | `DecimalFormat` with zero-check | Non-fisher formatting |

## Known Parity Traps

- **Trap A — Dual Formatting Scheme**: Non-fisher uses `format_f64()` (preserves trailing zeros). Fisher uses `format_rounded_value_to_print()` (strips trailing zeros). Same field → two different outputs depending on mode.
- **Trap B — hifreq Special Formatting**: In fisher mode, `hifreq` uses `format!("{:.4}")` NOT `getRoundedValueToPrint`. Preserves trailing zeros unlike other fisher fields.
- **Trap C — nm Clamping**: `nm` is clamped to 0 if negative, different formatting per fisher/non-fisher.
- **Trap D — Somatic msi Direct Formatting**: Somatic 61-col: msi uses `format!("{:.3}")` NOT `getRoundedValueToPrint`. Outputs "2.000" not "2".
- **Trap E — duprate Pattern Differs**: Simple non-fisher: "0.0" (1 decimal). Simple fisher: "0.00" (2 decimals). Somatic: "0.0" both.
- **Trap F — Null vs Placeholder Variants**: `empty()`: genotype="0", bias="0". `empty_null_variant()`: genotype="", bias="0;0".
- **Trap G — Somatic Constructor Parameter Ordering**: 4 params: begin_variant, end_variant, tumor_variant, normal_variant. Same object may be passed for multiple params.
- **Trap H — Amplicon Region Logic**: Non-null variant + good_variants → first good_variant's region. Null → "chr:pos-pos".
- **Trap I — Amplicon totalVariantsCount**: Uses `good_variants_count` (unique amplicons), NOT `good_variants.len()`.
- **Trap J — FisherExact odd_ratio**: Can return "Inf", integer string, or decimal. Must match Java formatting.
- **Trap K — Somatic Fisher P-value Direction**: Uses `min(p_value_less, p_value_greater)`.
- **Trap L — DecimalFormat HALF_EVEN Rounding**: Must use `round_half_even()` to match Java banker's rounding.
- **Trap M — Somatic Placeholder Detection**: Zero-count variant → special genotype/bias formatting.
- **Trap N — qratio Zero-Division**: No low-qual reads → multiply high_qual by 2 (not divide by 0.5).

## Divergences from Java

1. **No delimiter customization**: Rust hardcodes `\t`; Java `VariantPrinter` accepts a delimiter parameter.
2. **No VariantPrinter class**: Uses `Display` trait and standalone format functions instead.
3. **Fisher numerical precision**: Reimplements Apache Commons Math hypergeometric directly (no external dependency for the core math, uses `statrs` for CDF).
4. **Option<&Variant> vs null**: Distinct constructors (`empty()` / `empty_null_variant()`) replace Java null-handling.
5. **Debug computed externally**: `SomaticOutputVariant` debug flag set externally, not derived in constructor.

## Key Internal Functions

| Function | Purpose | Risk |
|---|---|---|
| `format_rounded_value_to_print()` | Fisher-mode formatter: whole-number detection + zero-stripping | HIGH |
| `format_f64()` | Non-fisher formatter: zero-check + trailing zeros preserved | HIGH |
| `format_fixed_f64()` | Direct format (debug output) | MEDIUM |
| `qratio_from_counts()` | High/low quality ratio with zero-division handling | HIGH |
| `is_somatic_placeholder_variant()` | Detects zero-count placeholder | MEDIUM |
| `format_somatic_genotype()` / `format_somatic_strand_bias()` | Placeholder-aware formatting | MEDIUM |
| `build_amplicon_debug()` | Debug string with tab/space delimiters | MEDIUM |
| `FisherExact::new()` | Hypergeometric setup | HIGH |
| `FisherExact::calculate_pvalue()` | Lower/upper/two-sided tails | HIGH |
| `FisherExact::mle()` | Odds ratio via root-finding | HIGH |
| `FisherExact::round_as_r()` | R-compatible rounding | HIGH |
| `zeroin_c()` | Brent's method root-finding | HIGH |
| `stirling_error()` / `deviance_part()` | Log-gamma approximation | HIGH |
| `log_hypergeometric_probability()` | PMF computation | HIGH |

## Cross-Module Dependencies

### Called By
- `vardict_pipeline.rs` — postprocessors create output variant structs and format to strings
- Tests

### Calls Into
- `to_vars_builder` — `Variant`, `VarType`, `StrandBiasFlag`, `var_type_string()`
- `scopedata/global_read_only_scope` — config (fisher, crispr, debug flags)
- `utils::round_half_even` — banker's rounding for Java `DecimalFormat` parity
- `statrs::distribution` — hypergeometric CDF (for `FisherExact`)
