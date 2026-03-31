---
name: repo-memory-distill
description: "Consolidate repository-scoped memory notes into thematic summaries. Use when: consolidate repo memory, distill memory files, merge memory notes, memory cleanup, reduce memory file count."
argument-hint: "Run when /memories/repo/ has more than ~50 files"
---

# Repo Memory Distill

## Purpose

Consolidate accumulated `/memories/repo/` files into thematic summary files for better scanability.

## When to Trigger

Trigger this workflow when `/memories/repo/` exceeds about 50 files. The current project has 59 files, so this skill is intended for the current scale of repository memory.

## Procedure

### Step 1: Inventory

- List all files in `/memories/repo/`.
- Read each file to understand its content.
- Note the date from the filename suffix, for example `_20260310`, when present.

### Step 2: Categorize by Theme

Adapt the categories to the actual content, but start from this mapping:

| Category | File | Content Types |
|----------|------|---------------|
| Parity Traps | `parity-traps.md` | Bugs that fooled diagnosis, misleading symptoms, false leads |
| Debugging Patterns | `debugging-patterns.md` | Effective debugging techniques, diagnostic approaches |
| Design Decisions | `design-decisions.md` | Architecture choices, why X over Y, API decisions |
| Known Limitations | `known-limitations.md` | Accepted divergences, platform-specific issues |
| Infrastructure | `infrastructure.md` | Build system, conda, test harness, CI notes |
| SV and Realigner | `sv-and-realigner.md` | Structural variant processor and realigner specifics |
| VecMap and Collections | `vecmap-and-collections.md` | Collection type decisions, ordering issues |
| Memory and Performance | `memory-and-performance.md` | Memory optimization, allocation, benchmarking |

Keep any files that do not fit these themes in their original form until a better category is clear.

### Step 3: Draft Summaries

For each category:

- Create a new summary file with a clear header.
- Group entries by sub-topic.
- Keep each entry to 1-3 bullet points max.
- Preserve the original date.
- Preserve unique insights and never discard information.
- Use this entry format:

```text
### <Topic> (YYYY-MM-DD)
- Key finding 1
- Key finding 2
- Related: <other entries or files if cross-referenced>
```

If a source note has no date suffix, keep the topic but mark the date conservatively, for example using the known file date if available elsewhere or `undated` if it cannot be recovered safely.

### Step 4: User Review

- Present the draft summaries to the user before deleting any originals.
- Show a mapping of which original files were merged into which summary.
- Get explicit approval before proceeding to deletion.

### Step 5: Replace

- Create the new summary files in `/memories/repo/`.
- Delete only the original files that were merged and approved for removal.
- Keep any files that do not fit into categories as-is.

## Quality Rules

- Never discard unique insights. If in doubt, keep the entry.
- Preserve original dates in each entry.
- Keep entries concise: 1-3 bullets per topic.
- Cross-reference related entries within summaries.
- If two entries contradict each other, keep both and note the contradiction.
- Reduce file count substantially, but keep the summaries lossless in information.

## Constraints

- Always get user approval before deleting any memory files.
- Never merge files across memory scopes: repo, session, and user memories must stay separate.
- Keep the total summary file count between 5 and 12, depending on content diversity.
- Each summary file should be scannable in under 30 seconds.