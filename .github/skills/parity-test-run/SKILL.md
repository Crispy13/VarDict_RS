---
name: parity-test-run
description: "Run parity tests for a specific subset of the config x chromosome matrix. Use when: run parity test, execute parity harness, test config, validate chromosome, parity sweep, smoke test."
argument-hint: "Specify scope like 'T1-01 chr20', 'smoke', 'dev', or 'release'."
---

# Parity Test Run Workflow

## Purpose

Run parity tests for a specific subset of the config x chromosome matrix.

## Constraints

- Do not modify source code while running tests.
- Always clear stale Rust cache before re-testing after a code change.
- Always include the harness exit code in the final report.
- Do not hardcode harness script names or fixed script paths. Discover the current harness in the workspace before executing anything.

## Procedure

### Step 1: Determine Scope

- **Single cell**: one config on one chromosome, for example `T1-01` on `chr20`
- **Preset**: `smoke` (9 cells), `dev` (100 cells), `release` (1100 cells)
- **Custom**: specific configs and/or chromosomes

Preset tiers:

| Preset | Config Scope | Chromosome Scope | Cells | Typical Use |
|--------|--------------|------------------|-------|-------------|
| smoke | 3 configs | 3 chromosomes | 9 | Quick sanity after small change |
| dev | 10 configs | 10 chromosomes | 100 | Overnight regression sweep |
| release | All 44 configs | All 25 chromosomes | 1100 | Full release validation |

Before running anything, restate the requested scope in `config x chromosome` form so the run is unambiguous.

### Step 2: Build

Build the Rust binary before testing.

- Use the `rust_build_env` conda environment.
- Use `cargo build --profile debug-release` for all development and verification builds.

Record which binary profile was built and which binary will be exercised by the harness.

### Step 3: Discover Harness

Search the workspace for parity test harness scripts instead of assuming a fixed script name.

Typical discovery targets:

- `tests/`
- `copilot-office/*/copilot-desk/scripts/`

Look for scripts that:

- Accept `--chr`, `--config-id`, `--preset`, `--opts` parameters
- Generate Java and Rust shard outputs and compare them
- Support `--rust-only`, `--no-stop`, `--parallel` flags

Discovery guidance:

- Search for executable shell or Python scripts mentioning `--config-id`, `--preset`, or `results.json`
- Confirm the script supports the requested scope before using it
- Prefer the harness that matches the current mission or parity workflow, but verify behavior from the script help text or source first

### Step 4: Execute

Use the discovered harness path in the commands below.

Common execution patterns:

**Post-fix verification** (most common):

```bash
# Clear stale Rust cache first
rm -rf tmp/na12878_parity/<label>/<chr>/rust/
rm -rf tmp/na12878_parity/<label>/<chr>/diff/

# Run with existing Java cache
<harness> --chr <N> --rust-only --no-stop --parallel 5
```

**Single config test**:

```bash
<harness> --config-id <ID> --chr <N> --no-stop --parallel 5
```

**Preset sweep**:

```bash
<harness> --preset <smoke|dev|release> --rust-only --no-stop
```

**Release validation**:

```bash
<harness> --preset release --rust-only --no-build --no-stop
```

Parallelism rules:

- Use `--parallel 5` for pileup (`-p`) configs because they are memory-intensive
- Use `--parallel 10` for non-pileup configs unless the discovered harness documents a stricter limit

Mode guidance:

- Use **cold** mode when Java artifacts must be regenerated
- Use **rust-only** mode when validating a Rust-only code change against existing Java cache

### Step 5: Collect Results

- Check the harness exit code. `0` means all requested cells passed.
- Look for `results.json` in `tmp/na12878_parity/<label>/<chr>/`
- Count `PASS`, `FAIL`, and `EMPTY` shards
- List failing shards and their diff files

If the harness spans multiple chromosomes or presets, aggregate results across all returned result directories before reporting.

### Step 6: Report

Use this report format:

```text
## Parity Test Report

**Scope**: {config(s)} x {chr(s)}
**Binary**: {debug-release | release}
**Mode**: {cold | rust-only}
**Harness Exit Code**: {code}

### Results
| Config | Chr | Shards | Pass | Fail | Empty |
|--------|-----|--------|------|------|-------|

### Failing Shards
| Config | Chr | Shard | Region | First Diff Line |
|--------|-----|-------|--------|-----------------|

### Next Steps
- Use `shard-diagnosis` skill on failing shards
- Use `mismatch-triage` skill to classify and prioritize
```

## Resource Constraints

- Use `--parallel 5` for pileup (`-p`) configs
- Use `--parallel 10` for non-pileup configs
- Monitor disk space because full-genome runs can consume 10-50 GB under `tmp/na12878_parity/`

## Known Pitfalls

- Stale Rust cache: always clear `rust/` and `diff/` directories after code changes
- Empty Java shards from OOM: delete 0-byte `.tsv` files before re-running
- Pileup memory pressure: keep parallelism capped at 5 for pileup mode