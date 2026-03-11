---
name: work_for_goal
description: Automates planning, execution, and reporting for the active mission autonomously.
---

<!-- 
This prompt does NOT require user confirmation between tasks, unlike `work_for_request`.
-->

You MUST perform the following plan-execute-report cycle:

## 1. Goal Planning
Run a custom agent named exactly 'Plan' as subagent to do the following:
1. Resolve `copilot-office/<mission-name>/` using `../instructions/project_context.instructions.md`.
	- Prefer the mission already established in the current session history.
	- If none exists, inspect the current request and closely related mission references to infer the mission safely.
	- Use `vscode_askQuestions` only if the mission still cannot be resolved.
2. Read the mission context from:
- `copilot-office/<mission-name>/copilot-project-plan.md`
- `copilot-office/<mission-name>/copilot-current-plan.md`
- `copilot-office/<mission-name>/copilot-desk/`
- `copilot-office/codebase/CODEBASE.md` when shared architecture matters
3. Ignore unrelated mission folders unless the user explicitly asks for cross-mission work.
4. Evaluate the current plan and progress in the active mission.
5. Update or generate `copilot-office/<mission-name>/copilot-current-plan.md` to maximize productivity while keeping tasks actionable and mission-scoped.

## 2. Autonomous Execution
You MUST implement the tasks defined in the active mission plan.
- When delegating to subagents, include the resolved mission name and tell them to stay within that mission folder and its related shared docs.
- Use a dedicated subagent for each distinct multi-step task (e.g., implementation, research, testing, debugging).
- Proceed sequentially through the active mission backlog WITHOUT stopping between tasks.
- Only pause execution if you encounter a critical blocker, require clarification, or hit a defined milestone requiring manual verification.
- Do not silently switch to another mission folder mid-run.
- Prefer direct tool calls in the main agent for trivial reads, single searches, and one-off commands.
- Use subagents when they reduce context pressure or isolate a bounded task; do not delegate every tool call mechanically.

## 3. Reporting
Once you reach a natural stopping point, Run a subagent to provide a comprehensive report detailing the resolved mission folder, whether it came from current session history or the inspection step, completed tasks, remaining backlog, and any uncovered issues.