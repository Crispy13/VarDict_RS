---
description: "Apply the multi-session workflow (goal hierarchy, stage cycling, persistent tracking) to any orchestrator agent. Adds session bootstrap, planning phases, auto-continuation, and work reporting without changing the agent's own implementation/review process."
---

# Multi-Session Workflow

You must follow the multi-session workflow defined below. This workflow governs how you manage goals, plan stages, track progress, and persist state across sessions. **Your own agent instructions define how you execute each stage** (your team, your implementation process, your review process). This prompt only defines the outer loop.

## How This Prompt Relates to Your Agent Instructions

This prompt and your agent's `.agent.md` serve different roles:

- **This prompt** = the **outer loop**. It controls *what* to work on next, *when* to start and stop a stage, and *how* to track progress across sessions.
- **Your agent instructions** = the **inner loop**. They control *how* you execute each stage — your team, your delegation workflow, your review process.

**The stage plan is the bridge.** When Phase 2 says "execute the stage," use your agent's own workflow to carry it out. The `copilot-stage-plan.md` defines the **scope boundary** — do not work outside it. When your workflow completes (your own review/test process approves the work), return to Phase 2's post-completion steps.

**Planning uses your available resources.** If your agent has a Planner subagent, use it for Phases 1a–1c. If not, plan directly using the tools you have.

**Drift Guard applies to both.** If your agent has a Drift Guard rule, re-read both your agent file and this prompt's rules after updating plan files.

---

## Goal Hierarchy

Organize work in three tiers using mission files:

| Tier | File | Scope |
|------|------|-------|
| **Project Goal** | `<mission>/copilot-project-plan.md` | The big-picture objective. Rarely changes. |
| **Active Goal** | `<mission>/copilot-active-plan.md` | The current milestone. Contains backlog, blockers, decisions, and context for resuming. |
| **Stage Goal** | `<mission>/copilot-stage-plan.md` | A single backlog item — the unit of work you execute using your own process. |

- `<mission>` = `copilot-office/missions/<mission-name>`
- **Project Desk** (`<mission>/copilot-desk/`) — decisions, specs, scratchpad notes.
- **Codebase Overview** (`copilot-office/codebase/CODEBASE.md`) — shared architecture and project structure.

---

## Workflow Phases

### Phase 0: Session Bootstrap

Run at the start of every session.

1. **Resolve mission folder.** On first request, ask the user for the mission name. Otherwise resolve from conversation history.
2. **Read context files** — only on first request or when you lack clear mission context. Skip if conversation already contains the information. Re-read a file only if it was modified since last read. Read in order:
   - `copilot-project-plan.md`
   - `copilot-active-plan.md`
   - `copilot-stage-plan.md` (if exists)
   - `copilot-office/codebase/CODEBASE.md`
3. **Determine entry point:**
   - Stage plan is in-progress → **Phase 2** (resume execution).
   - Stage plan is complete or missing → **Phase 1** (planning).
   - Active goal is complete → **Phase 3** (goal advancement).
   - No mission files → **Phase 1a** (project initialization).

### Phase 1: Planning

#### 1a: Project Initialization (first time only)

1. Ask the user about the project goal and scope.
2. Produce a project plan with milestones.
3. Write `copilot-project-plan.md`.
4. Create directory structure: `copilot-office/missions/<name>/copilot-desk/completed-stages/`.
5. Initialize `copilot-office/codebase/CODEBASE.md` if it doesn't exist.

#### 1b: Active Goal Selection

1. Ask the user to confirm: review `copilot-project-plan.md` and identify the next milestone, or set the active goal from the user's request.
2. Break the milestone into a working plan with backlog items, context, and acceptance criteria.
3. Write `copilot-active-plan.md`.

#### 1c: Stage Planning

1. Review `copilot-active-plan.md` to identify the next backlog item.
2. Produce a stage plan — a single focused sub-goal with: objective, context, files to modify/create, acceptance criteria, tests (if applicable).
3. Write `copilot-stage-plan.md`.

> **One backlog item = one stage.** Do not bundle multiple items.

### Phase 2: Stage Execution

**Execute the stage using your own agent's implementation and review process.** This prompt does not define how you implement — your agent instructions do. What this prompt requires:

- The stage plan (`copilot-stage-plan.md`) defines the scope.
- When execution completes (your process approves the work):

1. Mark `copilot-stage-plan.md` as `STATUS: COMPLETE`.
2. Update `copilot-active-plan.md`: mark backlog item done, note progress, update context.
3. Update `copilot-office/codebase/CODEBASE.md` if architectural changes were made.
4. Append a stage entry to `<mission>/copilot-work-report.md` (see Work Report format below).
5. Archive completed stage plan to `<mission>/copilot-desk/completed-stages/`.
6. Briefly notify the user (1-2 sentences), then:
   - More backlog items remain → **immediately** proceed to Phase 1c. Do NOT stop.
   - All items done → proceed to Phase 3.

> **Auto-continuation:** Keep looping through stages until the active goal is complete or the user interrupts.

### Phase 3: Goal Advancement

1. Mark `copilot-active-plan.md` as `STATUS: COMPLETE`.
2. Update `copilot-project-plan.md`: mark milestone done, note learnings.
3. Finalize the Work Report — add a summary section at the top.
4. Present milestone summary to the user.
5. Return to Phase 1b for the next active goal.

---

## File Formats

### copilot-stage-plan.md

```markdown
# Stage: {Title}

**Status:** IN-PROGRESS | COMPLETE
**Parent Goal:** {Active goal / backlog item reference}
**Created:** {Date}
**Last Updated:** {Date}

## Objective
{What this stage achieves}

## Context
{Relevant background, related files, constraints, prior decisions}

## Files to Modify/Create
- `path/to/file.ext` — {what needs to change}

## Acceptance Criteria
1. {Criterion 1}
2. {Criterion 2}

## Tests
{Commands to run, or "N/A"}

## Notes
{Blockers, decisions, gotchas}
```

### copilot-project-plan.md

```markdown
# Project: {Title}

**Mission:** {mission-name}
**Created:** {Date}
**Last Updated:** {Date}

## Vision
{High-level objective}

## Milestones
- [x] Milestone 1 — {Title} ✅
- [ ] Milestone 2 — {Title} ← ACTIVE
- [ ] Milestone 3 — {Title}

## Constraints & Principles
- {Key constraints}
```

### copilot-active-plan.md

```markdown
# Active Goal: {Title}

**Status:** IN-PROGRESS | COMPLETE
**Parent Project:** {Project reference}
**Created:** {Date}
**Last Updated:** {Date}

## Objective
{What this milestone achieves}

## Backlog
- [x] {Completed item}
- [ ] {Pending item} ← CURRENT STAGE
- [ ] {Future item}

## Blockers
- {Blocking issues}

## Decisions
- {Key decisions}

## Context
{What was just done, what's next, gotchas — enough to resume safely}
```

### copilot-work-report.md

```markdown
# Work Report: {Active Goal Title}

**Mission:** {mission-name}
**Started:** {Date}
**Last Updated:** {Date}
**Status:** IN-PROGRESS | COMPLETE

## Summary
{Added when goal completes}

---

## Stage 1: {Title}
**Date:** {Date}
**Backlog Item:** {Which item}

### What Was Done
- {Change — file, what, why}

### Files Changed
- `path/to/file.ext` — {description}

### Decisions Made
- {Trade-offs or decisions}

---

## Stage 2: {Title}
...
```

> **One report per active goal.** Archive previous to `copilot-desk/completed-reports/` when starting a new goal.

---

## Workflow Constraints

- **Persist relentlessly.** Every meaningful state change must be written to disk. A new session must resume cleanly from files alone.
- **Mission isolation.** Only access the resolved `<mission>` folder and `copilot-office/codebase/`. Never read other mission folders unless the user explicitly asks.
- **Respect the hierarchy.** One stage at a time. One backlog item per stage. Don't mix stage work with goal-level planning.
- **Auto-continue.** After completing a stage, immediately plan and execute the next one. Don't stop to report unless the active goal is done or the user interrupts.
