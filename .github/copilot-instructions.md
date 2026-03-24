# VarDictJava → Rust Parity Project

## Project Overview

This project is a Rust port of [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava), a variant discovery tool for next-generation sequencing (NGS) data. The primary goal is **100% output parity** with the original Java implementation while achieving idiomatic Rust and improved performance.

VarDictJava calls SNVs, MNVs, indels, complex variants, and structural variants from BAM files against a reference genome. It supports three modes: **Simple** (single-sample), **Somatic** (tumor/normal paired), and **Amplicon** (targeted sequencing).

## Architecture

### Java Pipeline (Reference)
```
SAMFileParser → RecordPreprocessor → CigarParser → VariationRealigner
    → StructuralVariantsProcessor → ToVarsBuilder → OutputVariant
```

### Key Modules and Parity Risk

| Module | Java LOC | Risk | Critical For |
|--------|----------|------|--------------|
| CigarParser | ~2,400 | HIGH | SNP/indel detection from CIGAR strings |
| VariationRealigner | ~1,200 | HIGH | Local realignment, position adjustment |
| StructuralVariantsProcessor | ~2,100 | HIGH | DEL/DUP/INV from discordant pairs |
| ToVarsBuilder | ~1,500 | MEDIUM | MSI, genotype, variant filtering |
| OutputVariant printers | ~400 | MEDIUM | Column formatting, float precision |
| FisherExact | ~200 | LOW | Strand bias statistics |
| Configuration/CLI | ~300 | LOW | Argument mapping |

## Parity Rules

1. **Output must be byte-identical** to Java for the same inputs — this is the non-negotiable standard
2. **Floating-point formatting**: Java uses `DecimalFormat("0.0000")`; Rust must match exactly, including trailing zeros and rounding behavior
3. **Collection ordering**: Use `IndexMap` instead of `HashMap` wherever Java uses `LinkedHashMap`
4. **Integer arithmetic**: Match Java's signed 32/64-bit overflow semantics using `wrapping_add`, `wrapping_mul` etc. where Java would silently overflow
5. **Null mapping**: Java `null` → Rust `Option::None`; every null-check branch in Java must have an equivalent `Option` match in Rust
6. **String handling**: Use `String`/`&str` with UTF-8; genomic data is ASCII-safe but validate at boundaries
7. **Tab-delimited output**: Columns must match exactly — count, order, and content

## Build and Test

```bash
# Build
cargo build --release

# Run tests
cargo test

# Run parity test against Java output
# (compare Rust output with reference Java output for test regions)
diff <(./target/release/vardict -G ref.fa -b test.bam -N sample regions.bed) expected_java_output.tsv

# Lint
cargo clippy -- -D warnings
cargo fmt --check
```

## Conventions

- Follow `rust.instructions.md` for general Rust style
- Every Rust function porting a Java method must reference the original Java class and method name in a doc comment
- When Java logic is intentionally complex or subtle, preserve the algorithmic structure even if it looks non-idiomatic — correctness over style
- Test each module independently against Java reference output before integration
- Use `#[cfg(test)]` modules co-located with implementation
- Thread model: `rayon` for data parallelism (replaces Java `CompletableFuture` + `ExecutorService`)
- BAM I/O: `rust-htslib` or `noodles` — document which is used and why

## Constraints
1. ~~Don't use `mimalloc`~~ — Fixed: The SIGSEGV was caused by `reallocate_c_vec` calling `System.dealloc` (libc free) on buffers allocated by mimalloc. Resolved by removing the C-malloc workaround after upgrading rust-htslib to v1.0.0c3. mimalloc is now the default allocator.