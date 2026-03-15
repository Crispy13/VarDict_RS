---
description: "Analyze a VarDictJava method for porting to Rust. Extracts algorithm logic, control flow, edge cases, null handling, and mutable state from Java source code."
# agent: "java-analyst"
argument-hint: "Specify the Java class and method, e.g. CigarParser.parseCigar"
---

Analyze the specified VarDictJava method for Rust porting.

Produce a complete analysis covering:

1. **Algorithm Logic**: Step-by-step description of what the method does
2. **Control Flow**: Every if/else, switch, loop, and early return with exact conditions
3. **Mutable State**: Every variable that changes during execution, with initial values and mutation points
4. **Null/Edge Cases**: Every null check, empty collection check, boundary condition
5. **Collection Ordering**: Whether any HashMap/LinkedHashMap iteration order affects output
6. **Float Formatting**: Any DecimalFormat usage and precision rules
7. **Dependencies**: External method calls and their return types
8. **Parity Warnings**: Specific concerns for Rust translation

Search for the Java source file in the workspace or reference the VarDictJava repository structure under `com.astrazeneca.vardict`.

Use the structured analysis format:
- Source location (file and line numbers)
- Pipeline stage position
- Parameter documentation
- Step-by-step algorithm with branch conditions
- Recommended Rust type mappings
