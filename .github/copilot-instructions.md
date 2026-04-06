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

## References

- Operational policies (build, env, test, sweep): see `.github/instructions/ops-policy.instructions.md`
- Parity rules (type mapping, float formatting, ordering): see `.github/instructions/rust-parity.instructions.md`
- Rust coding conventions: see `.github/instructions/rust.instructions.md`