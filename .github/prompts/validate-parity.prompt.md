---
description: "Validate output parity between VarDictJava and Rust port. Compares outputs byte-by-byte, identifies mismatches by column, and traces differences to their source."
# agent: "parity-tester"
argument-hint: "Specify what to validate: a module name, test region, or 'full' for end-to-end comparison"
---

Validate output parity for the specified scope.

Procedure:

1. **Identify test inputs**: Find or create appropriate BAM/BED test fixtures for the target module
2. **Generate reference output**: Run VarDictJava (if available) or locate saved reference output files
3. **Run Rust implementation**: Execute the Rust port with identical parameters
4. **Compare outputs**:
   - Line count comparison
   - Column-by-column diff for tab-delimited output
   - Float precision analysis for numeric fields
   - Order comparison for variant lines
5. **Report mismatches** in structured format:
   - Line number, column number, field name
   - Expected (Java) vs actual (Rust) values
   - Mismatch category (float formatting, missing variant, order, etc.)
6. **Trace root cause**: For the first mismatch, trace back through the code to identify the source

Output format:
- Overall result: PASS / FAIL with mismatch count
- Mismatch table with first 10 differences
- Root cause analysis for the primary issue
- Recommended next action (which method to re-analyze or fix)

For unit-level validation, create or update a Rust `#[test]` function that encodes the expected behavior.
