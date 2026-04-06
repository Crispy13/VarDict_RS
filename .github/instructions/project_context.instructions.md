---
description: Define locations and rules for project context, planning, and documentation files.
# applyTo: '**'
---

# Project Context Reference

The following files and directories provide project alignment and progress tracking. Resolve the active mission folder before making changes, and update relevant context files when milestones are reached. It is recommended to read these context files if you don't have context information.

## Definitions

- `<mission-folder>`: the resolved mission path placeholder. Set it to `copilot-office/missions/<mission-name>` before using any mission-specific files.


## 1. Project Goal (`<mission-folder>/copilot-project-plan.md`)
- **Purpose:** Contains the high-level objective and long-term plan for the active mission.


## 2. Active Plan (`<mission-folder>/copilot-active-plan.md`)
- **Purpose:** Contains the currently active goal, working backlog, blockers, unresolved decisions, and enough recent context to resume safely.


## 3. Stage Plan (`<mission-folder>/copilot-stage-plan.md`)
- **Purpose:** Contains the next actionable stage to execute for the active goal.


## 4. Codebase Overview (`copilot-office/codebase/CODEBASE.md`)
- **Purpose:** Contains shared architectural decisions, main components, and project structure across the workspace.


## 5. Project Desk (`<mission-folder>/copilot-desk/`)
- **Purpose:** A directory for storing mission-specific decisions, reference materials, technical specs, and scratchpad notes.
