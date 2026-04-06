---
name: tiered-config-test
description: "Graduated config testing workflow for VarDictJava-to-Rust parity. Use when: expanding config coverage, promoting between test tiers, config regression testing, pre-release validation, pairwise interaction testing."
argument-hint: "Specify tier like 'tier1', 'config-spread', or 'promote to core-wide'."
---

# Tiered Config Testing Workflow

## Harness Support

- The harness implements `smoke`, `dev`, `tier1`, `config-spread`, `core-wide`, `pairwise`, `full-gate`, and `release` presets.
- `BLOCKED_CONFIGS` and `--include-blocked` are implemented in the harness. The blocked list currently starts empty and should only grow when a config is intentionally quarantined.

## Purpose

Use a graduated testing pyramid to expand config coverage without defaulting to the full 44 x 25 release matrix on every change. Start with fast smoke validation, promote only after each gate passes, and reserve pairwise and full-release sweeps for broader interaction and ship-readiness checks.

## Testing Pyramid

| Tier | Preset Name | Configs | Chromosomes | Total Cells | Estimated Time | When to Run |
|------|-------------|---------|-------------|-------------|----------------|-------------|
| 0 | `smoke` | 3 (`T1-01`, `T1-03`, `T1-13`) | 3 (`20`, `22`, `MT`) | 9 | ~1 min | Every build, every code change |
| 1 | `tier1` | 14 (all `T1-*`) | 3 (`20`, `22`, `MT`) | 42 | ~5 min | After any parity fix |
| 2 | `config-spread` | 44 (all configs) | 3 (`20`, `22`, `MT`) | 132 | ~15 min | After CLI wiring changes, weekly |
| 3 | `core-wide` | 14 (all `T1-*`) | 25 (all real) | 350 | ~45 min | Pre-merge of major SV or pipeline changes |
| 4 | `pairwise` | 20 (`PW-000`..`PW-019`) | 10 (representative) | 200 | ~1-2 hr | Pre-release interaction testing |
| 5 | `release` | 44 (all configs) | 25 (all real) | 1100 | ~8-12 hr | Release candidate validation |

Note: the existing `dev` preset is 10 configs x 10 chromosomes = 100 cells, which places it between Testing Pyramid Tier `1` and Tier `2` in coverage and runtime.

Clarification: Testing Pyramid tiers `0` through `5` are execution stages defined by this skill. Harness `--tier N` values `1` through `4` are config-group filters inside the test matrix. They are not the same concept and should not be treated as interchangeable.

Composite preset:

| Preset | Chains | Behavior | Use |
|--------|--------|----------|-----|
| `full-gate` | `smoke` -> `tier1` -> `config-spread` -> `core-wide` | Early-exit on the first failing gate | Standard promotion path before pairwise or release |

Canonical rule: multi-config parity sweeps stop on the first failure, so preset runs should preserve that default behavior.

## Gate Criteria

- Tiers `0` through `3` require 100% `PASS` for all non-blocked configs. Any `FAIL` blocks promotion.
- Tier `4` (`pairwise`) requires at least 90% `PASS`. Treat new failures as backlog items unless the user explicitly elevates them to a release blocker.
- Tier `5` (`release`) requires 100% `PASS`. This is the release gate.

Gate process:

1. Run the requested tier preset.
2. Check the harness exit code.
3. If the run failed, diagnose the first failing shard with the `shard-diagnosis` skill.
4. Fix the root cause, rebuild, and re-run the tier.
5. Promote only after the tier meets its pass threshold.

## Promotion Procedure

1. Discover the active harness path in the workspace and confirm it supports the requested preset.
2. Build the Rust binary with `cargo build --profile debug-release` in the `rust_build_env` environment.
3. Run the current tier preset in `--rust-only` mode unless Java artifacts must be regenerated.
4. Check the harness exit code and inspect the results directory for failing shards.
5. If the tier passed, record the result in the active plan or mission notes and promote to the next tier.
6. If the tier failed, run the `shard-diagnosis` skill on the first failure and identify the root cause.
7. Implement the fix, rebuild the binary, and clear stale Rust and diff cache for the affected scope with `rm -rf tmp/na12878_parity/<label>/<chr>/rust/ tmp/na12878_parity/<label>/<chr>/diff/`.
8. Re-run Tier `0` (`smoke`) to catch regressions introduced by the fix.
9. Re-run the previously failing tier.
10. Repeat until the current tier passes, then continue promotion.

## Harness Commands

Discover the harness path before execution. Use `<harness>` as a placeholder for the discovered script path.

```bash
# Tier 0: Smoke
bash <harness> --preset smoke --rust-only

# Tier 1: All T1
bash <harness> --preset tier1 --rust-only

# Tier 2: Config spread (all 44, 3 chrs)
bash <harness> --preset config-spread --rust-only

# Tier 3: Core wide (T1, all chrs)
bash <harness> --preset core-wide --rust-only --parallel 5

# Tier 4: Pairwise (20 PW configs, 10 chrs)
bash <harness> --preset pairwise --rust-only --parallel 5

# Tier 5: Release
bash <harness> --preset release --rust-only --no-build

# Full gate (tiers 0 -> 1 -> 2 -> 3)
bash <harness> --preset full-gate --rust-only
```

### Manual Equivalents

The presets above are the standard entry points. For targeted debugging or custom scope control, the equivalent low-level harness flags are:

```bash
# Tier 1 interim: all Tier 1 configs on 20, 22, MT
bash <harness> --tier 1 --chr 20 --chr 22 --chr MT --rust-only

# Config-spread interim: all 44 configs on 20, 22, MT
bash <harness> --chr 20 --chr 22 --chr MT --rust-only

# Core-wide interim: all Tier 1 configs on all chromosomes
bash <harness> --tier 1 --all-chr --rust-only --parallel 5
```

- Pairwise: use `--preset pairwise` so the harness loads `tests/pairwise_configs.tsv` directly.
- Full-gate: use `--preset full-gate` for the standard chained promotion path, or run tiers `0` -> `1` -> `2` -> `3` manually when you need to inspect each gate separately.

## Pairwise Integration

- Pairwise configs live in `tests/pairwise_configs.tsv` and currently cover `PW-000` through `PW-019`.
- Generate or refresh the TSV with `tests/generate_pairwise_configs.py` when the option matrix changes.
- The TSV carries the exact `cli_flags` payload for each pairwise config.
- The harness should load `--preset pairwise` from the TSV at runtime rather than duplicating entries in shell arrays.
- Pairwise runs target option interactions that single-option configs and core-tier presets can miss.
- Fifteen of the twenty pairwise configs use pileup mode, so cap pairwise parallelism at `5`.

## Blocked Config Management

- `BLOCKED_CONFIGS` and `--include-blocked` are implemented in the harness.
- `BLOCKED_CONFIGS` currently starts empty, so all 44 configs are runnable by default.
- Add configs to `BLOCKED_CONFIGS` only when a failure is intentionally quarantined, and leave them visible in status output so the skipped scope stays explicit.

## Resource Constraints

| Setting | Value | Reason |
|---------|-------|--------|
| Pileup parallelism | 5 | Java can OOM above this on a 23 GB RAM machine |
| Non-pileup parallelism | 10 | Usually remains under 4 GB |
| Disk budget (full release) | 10-50 GB | Full-matrix artifacts and diffs accumulate quickly |
| Java cache | Preserve always | Deterministic and expensive to regenerate |
| Rust cache | Auto-invalidates on binary change | `.verified` markers track binary mtime |

## Known Pitfalls

- Stale Rust cache can invalidate conclusions. Clear `rust/` and `diff/` after code changes when re-testing the same scope.
- Empty Java shard outputs usually mean OOM or interrupted generation. Delete 0-byte `.tsv` files before re-running.
- Pileup-heavy presets can exceed memory budgets if parallelism drifts above `5`.
- Pairwise configs may expose unimplemented option interactions. Treat failures as real unless a config is explicitly quarantined in `BLOCKED_CONFIGS` after triage.
- `full-gate` stops on the first failing tier, so switch to the individual tier preset when you need finer debugging.