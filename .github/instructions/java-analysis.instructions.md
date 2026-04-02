---
description: "Use when analyzing VarDictJava source code, reviewing Java methods for porting, extracting algorithm logic from Java classes, or understanding Java control flow for Rust translation. Covers CigarParser, VariationRealigner, StructuralVariantsProcessor, ToVarsBuilder, and all VarDict modules."
applyTo: "**/*.java"
---

# VarDictJava Source Analysis Guidelines

## Analysis Objectives

When analyzing Java source from VarDictJava, extract the following for each method/class:

1. **Algorithm Logic**: Step-by-step description of what the code does, not just what it calls
2. **Control Flow**: Every `if/else`, `switch`, loop, and early `return` — these are parity-critical
3. **Mutable State**: Which variables are mutated, when, and how they affect downstream logic
4. **Edge Cases**: Null checks, boundary conditions, empty collections, zero-length strings
5. **Side Effects**: Modifications to parameters, shared state, or output streams
6. **Data Flow**: Input → transformation → output for each method

## VarDict-Specific Patterns to Watch

### Variant Description Strings
VarDictJava uses a dense string encoding for variants:
- `+SEQ` = insertion of SEQ
- `-N` = deletion of N bases
- `#SEQ` = complex variant
- `^` = intron
- `&` = reference match marker

Document every regex or string manipulation involving these patterns.

### Global State via `GlobalReadOnlyScope`
Java code accesses configuration through `GlobalReadOnlyScope.instance()`. Track every access — the Rust equivalent must pass configuration explicitly.

### HashMap Iteration
Java `LinkedHashMap` preserves insertion order. Document where iteration order matters for output.

### CIGAR Operations
Map each CIGAR operation (`M`, `I`, `D`, `N`, `S`, `H`, `=`, `X`) to the position arithmetic used. Off-by-one errors here destroy parity.

## Output Format

Structure your analysis as:

```
## Method: ClassName.methodName(params)

**Purpose**: One-line summary
**Called by**: Parent methods
**Calls**: Child methods

**Parameters**:
- param1: type — role in algorithm

**Algorithm**:
1. Step one...
2. Step two...
   - Branch A: condition → action
   - Branch B: else → action

**Mutable State Modified**:
- variable → how it changes

**Edge Cases**:
- When X is null/empty → behavior

**Parity Notes**:
- Specific concerns for Rust translation
```

## Priority Modules

Analyze in this order (by parity risk):
1. `modules/CigarParser.java` — core variant detection
2. `modules/VariationRealigner.java` — local realignment
3. `modules/StructuralVariantsProcessor.java` — SV detection
4. `modules/ToVarsBuilder.java` — variant building and filtering
5. `printers/*OutputVariant.java` — output formatting
6. `modes/SimpleMode.java`, `SomaticMode.java`, `AmpliconMode.java`

## Codebase Cache

A progressive Java codebase cache is maintained at `copilot-office/codebase/java/VarDictJava-CODEBASE.md`. Before doing a full source analysis of a Java method:

1. **Check the cache first**: Read the relevant module section in `VarDictJava-CODEBASE.md`. If the method has already been analyzed, verify the cached analysis against the actual Java source (spot-check key logic, don't re-read everything). Use the cached version if it's accurate.
2. **Analyze from source if missing**: If the method is not cached or the cache is incomplete for your needs, perform a full analysis from the Java source code.
3. **Update the cache after analysis**: Add your findings to the cache file following the Per-Module Template defined in the index. This is a secondary task — always return analysis results to the caller first.
