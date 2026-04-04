# VarDictJava Codebase Cache

> Progressive knowledge base built by `java-analyst` during source analysis.
> VarDictJava is **frozen at v1.8.3** — these analyses do not go stale.
> Each module has its own file to keep agent context small.

## Java Source Root

```
VarDictJava/src/main/java/com/astrazeneca/vardict/
```

## Pipeline Flow

```
Main → VarDictLauncher → Mode (Simple|Somatic|Amplicon)
  → SAMFileParser.parseSAM()
    → RecordPreprocessor.process()
    → CigarModifier.modifyCigar()
    → CigarParser.cigarParser()
    → VariationRealigner.realignVariation()
    → StructuralVariantsProcessor.processStructuralVariants()
    → ToVarsBuilder.toVars()
  → PostProcessModule
  → OutputVariant printer
```

## Module Index

| Module | Java File | LOC | Risk | Cache File | Status |
|--------|-----------|-----|------|------------|--------|
| CigarParser | `modules/CigarParser.java` | 2,662 | HIGH | [CigarParser.md](CigarParser.md) | complete |
| VariationRealigner | `modules/VariationRealigner.java` | 2,953 | HIGH | [VariationRealigner.md](VariationRealigner.md) | complete |
| StructuralVariantsProcessor | `modules/StructuralVariantsProcessor.java` | 2,307 | HIGH | [StructuralVariantsProcessor.md](StructuralVariantsProcessor.md) | complete |
| ToVarsBuilder | `modules/ToVarsBuilder.java` | 1,194 | MEDIUM | [ToVarsBuilder.md](ToVarsBuilder.md) | complete |
| CigarModifier | `modules/CigarModifier.java` | 787 | MEDIUM | [CigarModifier.md](CigarModifier.md) | complete |
| RecordPreprocessor | `modules/RecordPreprocessor.java` | 292 | LOW | [SAMFileParser.md](SAMFileParser.md) | complete |
| SAMFileParser | `modules/SAMFileParser.java` | ~342 | MEDIUM | [SAMFileParser.md](SAMFileParser.md) | complete |
| OutputVariant Printers | `printers/*.java` | 904 | MEDIUM | [OutputVariant.md](OutputVariant.md) | complete |
| Modes | `modes/*.java` | 691 | LOW | [Modes.md](Modes.md) | complete |
| Variations | `variations/*.java` | 1,125 | MEDIUM | [Variations.md](Variations.md) | complete |
| VariationMap | `collection/VariationMap.java` | 62 | LOW | [VariationMap.md](VariationMap.md) | complete |
| Configuration | `Configuration.java` | 374 | LOW | [Configuration.md](Configuration.md) | complete |
| Utils | `Utils.java` | 279 | LOW | [Utils.md](Utils.md) | complete |
| ScopeData | `data/scopedata/*.java` + `data/*.java` | ~950 | LOW | [ScopeData.md](ScopeData.md) | complete |

## Per-Module File Template

When creating a module file for the first time, use this structure:

```markdown
# ModuleName

**Source**: `path/relative/to/vardict/`
**LOC**: N
**Rust counterpart**: `src/mods/module_name.rs`
**Status**: partial | complete

## Overview
One-paragraph summary of the module's role in the pipeline.

## Method Inventory
| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| methodName() | L100-L200 | yes/no | one-liner |

## Method Analyses
(Full step-by-step analyses following the java-analyst output format)

## Cross-Module Dependencies
- Calls: which modules this depends on
- Called by: which modules invoke this

## Known Parity Traps
- Specific Java behaviors that cause Rust translation bugs
```
