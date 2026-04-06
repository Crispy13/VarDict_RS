---
name: codebase-doc-manage
description: "Manage the Java and Rust codebase documentation cache. Use when: orient on a module before analysis, update codebase docs after a fix, audit cache completeness, check module structure, reduce cold-start reading time for stateless subagents."
argument-hint: "Module name to orient on, or 'audit' to check cache completeness"
---

# Codebase Documentation Management

## Purpose

Maintain the progressive codebase caches so **stateless subagents** can orient on any module in seconds rather than re-reading thousands of source lines per invocation.

The caches answer three questions for any parity issue:
1. **Where to read** — which files, which line ranges, which methods
2. **How the program is structured** — pipeline flow, module boundaries, data flow
3. **What traps exist** — parity-critical behaviors already discovered

## Cache Locations

```
copilot-office/codebase/
├── java/
│   ├── VarDictJava-CODEBASE.md    ← Index (pipeline flow, 14-module table, template)
│   └── {ModuleName}.md            ← Per-module doc, created on first analysis
└── rust/
    ├── VarDict-rs-CODEBASE.md     ← Index (pipeline flow, 14-module table, template)
    └── {module_name}.md           ← Per-module doc, created on first analysis
```

## Phase 1: Orient (Read)

Use this at the start of any task to avoid re-reading source from scratch.

1. **Read the index** for the relevant language:
   - Java: `copilot-office/codebase/java/VarDictJava-CODEBASE.md`
   - Rust: `copilot-office/codebase/rust/VarDict-rs-CODEBASE.md`
   Read both if the task involves Java→Rust correspondence.

2. **Find the target module row** in the index table. Note its Cache File link and Status.

3. **If Status is `partial` or `complete`**: Read the module doc file. It contains:
   - Overview and pipeline role
   - Method inventory (which methods are analyzed)
   - Known parity traps (directly relevant to your task)
   - Java↔Rust correspondence (for cross-language work)

4. **If Status is `not started`**: No module doc exists yet. Proceed to source analysis.
   After analysis, run Phase 2.

5. **Use what you find** — don't re-analyze methods already documented unless you need to verify a specific parity trap.

## Phase 2: Write / Update

Use this after completing analysis or implementation to preserve what you learned.

**Creating a new module doc** (Status was `not started`):
1. Read the **Per-Module File Template** from the relevant index file.
2. Create the module doc file at:
   - Java: `copilot-office/codebase/java/{ModuleName}.md`
   - Rust: `copilot-office/codebase/rust/{module_name}.md`
3. Populate the sections using what you analyzed. Minimum required:
   - Overview (one paragraph)
   - Method Inventory table (fill in what you analyzed; mark others as `no`)
   - Known Parity Traps (even if empty — write "None found yet")
4. Update the index table: change Status from `not started` → `partial` or `complete`.

**Updating an existing module doc** (Status was `partial` or `complete`):
1. Read the existing module doc.
2. Touch ONLY the sections changed by your current work:
   - New method analyzed? Add row to Method Inventory (mark `yes`). Add method analysis.
   - New parity trap found? Append to Known Parity Traps. Never delete existing traps.
   - Rust divergence discovered? Update Divergences section.
   - Java↔Rust correspondence changed? Update correspondence table.
3. Do NOT rewrite existing correct content.
4. If the module is now fully documented, update Status in the index to `complete`.

**Content rules:**
- **Java docs** — detailed is good (source is frozen, won't stale). Include step-by-step method analyses.
- **Rust docs** — architecture-level only (code evolves). Record Java correspondence, parity traps, and divergences. Skip method-level Rust detail.
- **Both** — always include Cross-Module Dependencies and Known Parity Traps.
- **Keep it concise** — this is a reference, not a source dump. One paragraph per method overview; full detail only for parity-critical logic.

## Phase 3: Audit

Run when asked to check cache health or before a large parity sweep.

1. Read both index files.
2. For each module row, check Status:
   - `not started` → **Gap**. Flag if Risk is HIGH or MEDIUM.
   - `partial` → Read the module doc. Count `no` entries in Method Inventory. Flag HIGH-risk modules with >50% unanalyzed methods.
   - `complete` → Spot-check: does the doc reflect current source structure? (Especially Rust docs.)
3. Output an audit report:

```
## Codebase Cache Audit — {date}

### Java Cache
| Module | Risk | Status | Gap? |
|--------|------|--------|------|
| CigarParser | HIGH | not started | YES |

### Rust Cache
| Module | Risk | Status | Gap? |
|--------|------|--------|------|
| cigar_parser | HIGH | not started | YES |

### Priority Recommendations
1. {module} — {risk} risk, {reason to document first}
```

## Documentation Gate Protocol

Use this after receiving a report from `java-analyst` or `rust-implementer`.

1. Save the report to a session file: `/memories/session/{agent}-{module}-report.md`.
2. Dispatch `codebase-librarian` with `report_path`, `module`, `language`, and `mode: update`.
3. Verify the librarian response contains a `Cache Update:` footer before you treat the doc task as complete.

| Footer Value | Action |
|--------------|--------|
| `Cache Update: wrote ...` | Proceed — cache populated |
| `Cache Update: updated ...` | Proceed — cache extended |
| `Cache Update: no actionable content ...` | Acceptable — log and proceed |
| Footer absent or error | Re-dispatch the librarian once. If it fails twice, log the gap and continue; do not block parity work on doc failures. |

Before switching to a new module, run a module-transition audit by dispatching `codebase-librarian` in `audit` mode for modules touched in the current session.

## Constraints

- Only edit files under `copilot-office/codebase/java/` and `copilot-office/codebase/rust/`.
- Never modify source code files (`.rs`, `.java`).
- Always update the index Status when creating or substantially updating a module doc.
