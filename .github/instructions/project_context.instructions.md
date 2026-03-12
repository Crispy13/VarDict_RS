---
description: Define locations and rules for project context, planning, and documentation files.
# applyTo: '**'
---

# Project Context Reference

You must utilize the following files and directories to maintain project alignment and track progress. Resolve the active mission folder before making changes, and update relevant context files when milestones are reached.

## 0. Mission Folder Resolution (`copilot-office/<mission-name>/`)
- **Purpose:** Defines the active mission workspace under `copilot-office/`.
- **Agent Rule:** Resolve `<mission-name>` in this order: use the mission already established in the current session history if one exists; otherwise, inspect the current user request and closely related mission references to infer it; otherwise, ask the user to choose the mission with `vscode_askQuestions`.
- **Agent Rule:** Never silently invent `<mission-name>` or keep the literal placeholder in real work.
- **Agent Rule:** Once a mission is resolved in the current session, keep using that mission unless the user explicitly changes it or asks for cross-mission work.
- **Agent Rule:** The inspection step should stay narrow: prefer explicit mission names in the current request, recent chat context, or directly related mission documents. Do not infer from unrelated open files or weak codebase similarity.
- **Agent Rule:** Ask the user only when the current session does not already establish the mission and the inspection step still cannot resolve it safely.
- **Agent Rule:** Ignore unrelated mission folders unless the user explicitly asks for cross-mission work.

## 1. Project Goal (`copilot-office/<mission-name>/copilot-project-plan.md`)
- **Purpose:** Defines the high-level objective and long-term plan for the active mission.
- **Agent Rule:** Read this first to understand the active mission. Propose updates only when the mission direction or scope changes.
- **Agent Rule:** If mission-specific files conflict with generic notes, follow the mission-specific files and surface the conflict.

## 2. Active Plan (`copilot-office/<mission-name>/copilot-active-plan.md`)
- **Purpose:** Tracks the currently active goal for the mission, including the working backlog, blockers, unresolved decisions, and enough recent context to resume safely.
- **Agent Rule:** Read this to understand what goal is currently active and what still needs to be achieved.
- **Agent Rule:** Update this file when the active goal changes materially or when meaningful progress changes the working state.
- **Agent Rule:** Keep this file operational and current. Compress stale history and move durable design notes into the mission desk rather than turning the active plan into an archive.

## 3. Stage Plan (`copilot-office/<mission-name>/copilot-stage-plan.md`)
- **Purpose:** Defines the next actionable stage to execute for the active goal.
- **Agent Rule:** Read this file to know what to execute now.
- **Agent Rule:** Regenerate or rewrite this file when a stage is completed, replaced, or re-scoped.
- **Agent Rule:** Keep this file short, execution-oriented, and specific enough to drive immediate work.

## 4. Codebase Overview (`copilot-office/codebase/CODEBASE.md`)
- **Purpose:** Explains shared architectural decisions, main components, and project structure across the workspace.
- **Agent Rule:** Use this file to understand shared architecture and conventions.
- **Agent Rule:** Do not let this global file override mission-specific requirements in the active mission folder.
- **Agent Rule:** Keep this file limited to workspace-wide architecture, common modules, shared infrastructure, and conventions reused by multiple missions.
- **Agent Rule:** Do not store mission-specific backlog, acceptance criteria, temporary status, or mission-only decisions here. Those belong in the active mission folder.

## 5. Project Desk (`copilot-office/<mission-name>/copilot-desk/`)
- **Purpose:** A directory for storing mission-specific decisions, reference materials, technical specs, and scratchpad notes.
- **Agent Rule:** Save important design choices, data schemas, and reference snippets here as separate Markdown files to retain knowledge across sessions.