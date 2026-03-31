---
name: parity-workflow
description: "Operational guide for the VarDictJava-to-Rust parity infrastructure. Use when: parity testing workflow, run parity harness, cache management, golden archive, shard comparison, test infrastructure, parity verification, config matrix."
argument-hint: "Describe what aspect of the parity workflow you need help with"
---

# Parity Workflow

## Purpose

Operational guide for the VarDictJava -> Rust parity test infrastructure.

This skill covers how to run and verify parity infrastructure safely: harness discovery, cache semantics, shard execution, presets, archive usage, and common failure modes.

## Scope Boundary

- Use this skill for runtime operations and parity harness behavior.
- Do not use this skill for the full analyze -> implement -> test -> review loop.
- Hand code-fix workflows to `parity-check` after the failing shard or failing field is identified.

## Infrastructure Overview

Core script roles:

- `shard runner`: runs one config on one chromosome, generates Java and Rust shard outputs, and compares them.
- `option matrix driver`: runs many config x chromosome cells and aggregates results.
- `status reporter`: summarizes cached pass, fail, and empty states without re-running shards.
- `cleanup`: removes stale Rust and diff artifacts and optionally prunes broader caches.
- `golden archive`: compresses and restores reusable Java cache for faster re-runs.

Script locations are not fixed forever. Discover current entrypoints under:

- `tests/`
- `copilot-office/*/copilot-desk/scripts/`

Do not hardcode a single path when executing the harness in automation or instructions.

## Cache Directory Structure

```text
tmp/na12878_parity/
├── <opts-label>/
│   └── <chr>/
│       ├── java/shard_NNN.tsv    # Deterministic, reusable
│       ├── rust/shard_NNN.tsv    # Auto-invalidated on binary change
│       ├── diff/shard_NNN.diff   # Unified diff
│       ├── diff/shard_NNN.meta   # Structured metadata
│       └── diff/shard_NNN.status # PASS|FAIL|EMPTY
```

Cache rules:

- Java outputs are preserved on cleanup because they are deterministic and expensive to regenerate.
- Rust outputs are invalidated and regenerated when the Rust binary is newer than the cached shard output.
- `.verified` markers persist until the Rust binary changes.
- `--rust-only` skips Java generation and requires an existing valid Java cache.
- `results.json` may exist alongside shard directories and is the main machine-readable summary for a cell.

Operational guidance:

- After a Rust code change, remove `rust/` and `diff/` before trusting a re-run.
- Preserve `java/` unless the Java cache is corrupt, incomplete, or intentionally being rebuilt.
- Treat empty Java shard files as broken cache, not successful prior work.

## Shard Runner Execution Model

The shard runner is a 4-phase workflow:

1. Java generation
2. Rust generation
3. Comparison
4. Cleanup

Execution guarantees:

- Java and Rust do not run simultaneously. This is deliberate and prevents memory exhaustion.
- Comparison is byte-for-byte and writes diff, meta, and status artifacts.
- Cleanup usually removes passing Rust and diff outputs to save disk while retaining Java cache.

Shard semantics:

- Shards are 1 MB genomic regions.
- A shard boundary stays fixed even when region extension is enabled.
- Region extension via `-x` expands internal fetch windows, not the shard coordinate range.
- Boundary reads can still change visible output inside a shard because the fetch window is wider than the shard itself.

## Config Tiers And Presets

Quick tier reference:

- `T1` (14): core modes such as pileup, nosv, freq, realign, mapq, fisher, debug.
- `T2` (16): single-parameter variations such as qual, vext, mismatch, filter, readpos.
- `T3` (8): multi-parameter combinations.
- `T4` (6): region and size parameters.

Preset scopes:

| Preset | Cells | Typical Use |
|--------|-------|-------------|
| `smoke` | 9 | Fast sanity after a small change |
| `dev` | 100 | Broader development or overnight sweep |
| `release` | 1100 | Full validation |

Preset behavior:

- Presets define both config scope and chromosome scope.
- Per-dimension flags such as `--chr` can narrow a preset after selection.
- Use presets for repeatable verification, not ad hoc coverage.

## Known Pitfalls

These are the critical failure modes that waste the most time.

### 1. Stale Cache Trap

After a Rust code change, cached Rust shards can survive and make a fix appear ineffective.

Required action:

```bash
rm -rf tmp/na12878_parity/<label>/<chr>/rust/
rm -rf tmp/na12878_parity/<label>/<chr>/diff/
```

Always do this before re-testing a fix.

### 2. Region Extension Amplifies Bugs

`-x 150` widens the internal fetch window and pulls in boundary reads. This can surface bugs that do not appear in lighter configs.

Implication:

- A config that passes without extension can still fail with extension.
- Heavy configs are often better bug finders than default configs.

### 3. Pileup Memory Exhaustion

Pileup mode is materially heavier than non-pileup runs.

Rule:

- Use `--parallel 5` for pileup configs.
- Do not use `--parallel 10` for pileup on typical development hardware.

### 4. Empty Java Cache From OOM

If Java fails mid-shard, a 0-byte `.tsv` can remain in cache and poison future runs.

Required action:

```bash
find tmp/na12878_parity/<label>/<chr>/java/ -name 'shard_*.tsv' -empty -delete
```

Then re-run the affected cell.

### 5. Direct `-R` Run Requires Deterministic Java

When reproducing a shard outside the harness, Java must run with `-th 1`.

Rule:

- Use `-R "chr:start-end"` for both implementations.
- Add `-th 1` to Java for deterministic output.

### 6. Shard Boundary Effects

The same variant can appear with different apparent coverage or neighboring context near shard edges.

Interpretation rule:

- Check whether a mismatch is a true logic difference or a boundary-sensitive manifestation of fetch extension and shard partitioning.

## Verification Checklist

Before declaring a parity fix complete:

- [ ] `cargo test` reports 0 failures.
- [ ] The target config passes on the target chromosome.
- [ ] A default non-pileup regression check still passes.
- [ ] Stale Rust cache was cleared before the validation run.
- [ ] Debug probes or temporary diagnostics were removed from source.
- [ ] Plan files were updated.

Interpretation:

- Passing the target cell alone is not enough.
- A cache-contaminated pass is not trustworthy.
- Verification should be stated in `config x chromosome` terms so the claim is auditable.

## Golden Archive Workflow

Purpose:

- Share reusable Java cache.
- Restore known-good Java outputs without paying regeneration cost again.
- Speed up Rust-only re-validation after code changes.

Workflow:

1. Create an archive from Java cache as `.tar.zst`.
2. Restore the archive on another machine or after cleanup.
3. Re-run parity in `--rust-only` mode.

Operational notes:

- Pileup is excluded by default because it is much larger.
- Non-pileup Java cache is small enough to archive routinely.
- Archive restore is useful after disk cleanup or when sharing a stable baseline across sessions.

Typical size expectations:

- Non-pileup archive: about 200 MB compressed.
- Pileup archive: usually omitted by default because it can reach tens of GB.

## Resource Constraints

| Setting | Value | Reason |
|---------|-------|--------|
| `--parallel 5` | Pileup configs | Java can exceed 10 GB RSS with 10 workers |
| `--parallel 10` | Non-pileup configs | Usually stays under about 4 GB |
| Timeout | 600s per shard | Default harness timeout |
| Disk | 10-50 GB | Full-genome parity runs, especially with pileup |

Operational consequences:

- Memory limits are usually driven by Java generation, not Rust comparison.
- Disk pressure accumulates under `tmp/na12878_parity/` and should be monitored during broad sweeps.
- Cached Rust-only re-runs are far cheaper than cold runs, but still require disk headroom for regenerated Rust and diff outputs.

## Execution Guidance

Preferred modes:

- Use cold runs when Java cache is missing or intentionally being rebuilt.
- Use `--rust-only` when validating a Rust-only code change against existing Java cache.
- Use release binaries for longer verification sweeps when runtime matters.

Safety rules:

- Build the intended Rust binary before a large run.
- Do not assume the harness will choose the correct binary profile unless its flags explicitly say so.
- If cached artifacts and binary timestamps are suspicious, clear Rust cache explicitly instead of guessing.

## Related Skills

- `parity-check`: full analyze -> implement -> test -> review workflow after the operational failure is localized.
- `shard-diagnosis`: investigate a specific failing shard and identify the first divergent line or field.
- `mismatch-triage`: classify and prioritize many mismatches.
- `parity-test-run`: execute parity tests for a specific scope such as a cell, preset, or sweep.

## Handoff Guidance

Choose the next skill based on the current state:

- Harness or cache question -> stay in `parity-workflow`.
- Need to run a specific test scope -> use `parity-test-run`.
- Need to inspect one failing shard -> use `shard-diagnosis`.
- Need to prioritize a larger batch of mismatches -> use `mismatch-triage`.
- Need to fix the underlying implementation -> use `parity-check`.