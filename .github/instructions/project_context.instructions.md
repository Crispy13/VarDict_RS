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

## 2. Current Working Plan (`copilot-office/<mission-name>/copilot-current-plan.md`)
- **Purpose:** Tracks immediate next steps, current tasks, and recent progress for the active mission.
- **Agent Rule:** Read this to know what to work on next. Update this file continuously as tasks are completed or new actionable steps are identified.
- **Agent Rule:** Keep this file operational and current. Move durable design notes into the mission desk rather than turning the current plan into an archive.

## 3. Codebase Overview (`copilot-office/codebase/CODEBASE.md`)
- **Purpose:** Explains shared architectural decisions, main components, and project structure across the workspace.
- **Agent Rule:** Use this file to understand shared architecture and conventions.
- **Agent Rule:** Do not let this global file override mission-specific requirements in the active mission folder.
- **Agent Rule:** Keep this file limited to workspace-wide architecture, common modules, shared infrastructure, and conventions reused by multiple missions.
- **Agent Rule:** Do not store mission-specific backlog, acceptance criteria, temporary status, or mission-only decisions here. Those belong in the active mission folder.

## 4. Project Desk (`copilot-office/<mission-name>/copilot-desk/`)
- **Purpose:** A directory for storing mission-specific decisions, reference materials, technical specs, and scratchpad notes.
- **Agent Rule:** Save important design choices, data schemas, and reference snippets here as separate Markdown files to retain knowledge across sessions.