---
name: parity-check
description: "Full parity check workflow for VarDictJava-to-Rust port. Use when running a complete analyze-implement-test-review cycle for a Java module, fixing parity mismatches end-to-end, or verifying output parity for VarDict variant calling."
argument-hint: "Specify the Java module or method to check, e.g. 'CigarParser' or 'CigarParser.parseCigar'"
---

# Parity Check Workflow

End-to-end workflow for achieving and verifying 100% output parity between a VarDictJava module and its Rust port.

## When to Use

- Porting a new Java module/method to Rust
- Fixing a known parity mismatch
- Verifying a module after refactoring
- Full audit of a Rust module against its Java original

## Prerequisites

- VarDictJava source code accessible (local or via web)
- Rust port codebase with `Cargo.toml`
- Test BAM/BED fixtures (or ability to create minimal fixtures)
- Reference Java output for comparison (or ability to run Java version)

## New Module Track

### Phase 1: Analysis

1. **Identify the target**: Determine which Java class/method needs parity work
2. **Delegate to java-analyst**: Request structured analysis of the Java method covering:
   - Algorithm logic (step-by-step)
   - Control flow (every branch)
   - Mutable state tracking
   - Null/edge cases
   - Collection ordering dependencies
   - Float formatting patterns
3. **Review the analysis**: Verify completeness — every branch, every null check, every side effect must be documented
4. **Identify parity risk areas**: Mark HIGH/MEDIUM/LOW risk for each aspect

### Phase 2: Implementation

5. **Read existing Rust code**: Check if a partial implementation exists
6. **Delegate to rust-implementer**: Provide the java-analyst output and request Rust implementation following:
   - Type mapping from `rust-parity.instructions.md`
   - Traceability comments linking to Java source
   - Every branch from the analysis implemented
   - `cargo check` verification
7. **Verify compilation**: Ensure the code compiles cleanly
7.5. **Pre-test structural review**: Before running the parity sweep, verify:
    - Every Java branch from the analysis has a corresponding Rust branch (no TODO/unimplemented stubs)
    - All output-affecting maps use `IndexMap` (not `HashMap`) where Java uses `LinkedHashMap`
    - Float formatting helpers match Java `DecimalFormat` precision rules
    If any check fails, loop back to Phase 2 before spending time on the expensive parity sweep.

### Phase 3: Testing

8. **Prepare test cases**: Identify or create test inputs that exercise:
   - Normal case (typical variant calling scenario)
   - Edge cases from the analysis (empty input, null equivalents, boundary values)
   - Numeric precision cases (values that test float rounding)
9. **Delegate to parity-tester**: Run comparison between Java reference output and Rust output
10. **Analyze results**:
    - PASS → proceed to Phase 4
    - FAIL → examine mismatch report, loop back to Phase 1 for the failing area

### Phase 4: Review

11. **Delegate to code-reviewer**: Request quality review covering:
    - Parity correctness checklist
    - Idiomatic Rust compliance
    - Performance impact assessment (using `change-impact-review` skill)
    - Extensibility evaluation
12. **Verify Performance Verdict**: Confirm the code-reviewer produced a Performance Verdict:
    - `PERF_SAFE` → proceed
    - `PERF_RISK` → document and notify user
    - `PERF_REGRESSION` → STOP, escalate before proceeding
13. **Address review findings**:
    - Blocking issues → loop back to Phase 2
    - Non-blocking suggestions → create follow-up tasks
14. **Mark module as verified**: Update tracking with parity status

## Bug Fix Track

1. **Reproduce**: Capture the specific output difference (expected vs actual).
2. **Trace**: Delegate to `java-analyst` to trace the mismatching field or behavior back to the source logic.
3. **Write failing test**: Extract the fixture from Java output (never from Rust). Add a `#[test]` in the appropriate fixture test file that is `#[ignore]`d and fails before the fix. See the **Find → Fix → Test** rule in `rust-parity.instructions.md` for naming conventions and fixture requirements. This test is the acceptance gate for the fix.
4. **Fix**: Delegate to `rust-implementer` with the exact Java logic, mismatch details, and the test from step 3 as the acceptance gate.
5. **Doc Gate**: Save the producer reports, then follow `codebase-doc-manage` Phase 2 and `Documentation Gate Protocol` for cache updates.
6. **Pre-Fix Review**: Delegate to `code-reviewer` for Parity Correctness only — verify the test from step 3 exists, was failing before the fix, and now passes. If changes are requested, loop back to step 4.
7. **Validate**: Delegate to `parity-tester` to confirm the fix on the affected case or shard.
8. **Post-Fix Review**: Delegate to `code-reviewer` for the full review (all four sections, including the binding Performance Verdict).
9. **Performance Gate**: Act on the verdict: `PERF_SAFE` → proceed, `PERF_RISK` → log and notify, `PERF_REGRESSION` → block and escalate.

## Iteration Protocol

If parity test fails:
1. Take the FIRST mismatch from the test report
2. Trace it to the specific Java logic branch
3. Re-analyze ONLY that specific area (not the whole method)
4. Fix the specific Rust code path
5. Re-run parity test
6. Repeat until all mismatches are resolved

## Module Priority Order

Process VarDictJava modules in this order (highest parity risk first):

| Priority | Module | Java LOC | Key Challenge |
|----------|--------|----------|---------------|
| P0 | CigarParser | ~2,400 | Complex mutable state, CIGAR operations |
| P0 | VariationRealigner | ~1,200 | Float arithmetic, position adjustments |
| P0 | StructuralVariantsProcessor | ~2,100 | Discordant pairs, split reads |
| P1 | ToVarsBuilder | ~1,500 | MSI detection, variant filtering |
| P1 | *OutputVariant printers | ~400 | Column formatting, float precision |
| P2 | SimpleMode / SomaticMode / AmpliconMode | ~800 | Pipeline orchestration |
| P3 | FisherExact | ~200 | Mathematical (well-defined) |
| P3 | Configuration / CmdParser | ~300 | Straightforward mapping |

## Success Criteria

A module passes parity check when:
- [ ] Rust output is byte-identical to Java output for all test inputs
- [ ] All edge cases from the analysis are tested
- [ ] Code review has no blocking issues
- [ ] Performance Verdict is `PERF_SAFE` or `PERF_RISK` (with documented justification)
- [ ] `cargo test -- --include-ignored` passes for the module
- [ ] `cargo clippy` has no warnings for the module
