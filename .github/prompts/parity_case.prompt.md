---
name: parity_case
description: "Investigate a parity manifest row, case file, BAM/reference pair, or tag and return the most likely failing stage plus the narrowest validation path."
argument-hint: "Manifest row, case file, tag, or BAM/reference pair"
agent: "parity-investigator"
---

Investigate the requested parity case using [tests/parity_case_manifest.csv](../../tests/parity_case_manifest.csv), [tests/integration_test.rs](../../tests/integration_test.rs), and the relevant pipeline modules under [src/mods](../../src/mods).

Your goal is diagnosis, not broad implementation.

Return:
- The exact manifest row or filtered case set you matched.
- The mode, options, reference, BAM, and tags that matter for this case.
- The most likely failing pipeline stage or boundary.
- The strongest code or test evidence supporting that diagnosis.
- The narrowest next validation step to confirm the hypothesis.
- Any open questions or missing context that block a confident fix.

If multiple cases match, narrow to the smallest meaningful set and explain why.