---
name: work_for_goal
description: Automates planning, execution, and reporting for the active mission autonomously.
---

<!-- 
This prompt does NOT require user confirmation between tasks, unlike `work_for_request`.
-->

# Workflow
You MUST perform the following workflow:

## Definitions
- `Mission goal`: the higher-level objective described in `copilot-office/<mission-name>/copilot-project-plan.md`.
- `Active plan`: `copilot-office/<mission-name>/copilot-active-plan.md`, which tracks the currently active goal, the working state needed to continue it, and the high-level slices or expected stages for that goal.
- `Stage plan`: `copilot-office/<mission-name>/copilot-stage-plan.md`, which defines the next actionable stage to execute.
- `Active goal`: the concrete goal you are currently trying to achieve in this run. It may be a user-provided subgoal or, if no subgoal is active, the mission goal itself.
- `Stage objective`: the specific execution slice of the active goal described in the current stage plan.

## 0. Resolve The Active Goal
- First resolve `copilot-office/<mission-name>/` using `../instructions/project_context.instructions.md`.
- Then determine the `active goal` using this priority:
1. If the user gives an explicit request in the current prompt, that request is the active goal.
2. Otherwise, if session history or `copilot-office/<mission-name>/copilot-active-plan.md` already shows an unresolved active goal, continue that active goal.
3. Otherwise, use the mission goal from `copilot-office/<mission-name>/copilot-project-plan.md` as the active goal.
- Do not confuse the mission goal with the active goal. The mission goal is the larger objective; the active goal is what you are currently pursuing.
- If the active goal is vague or cannot be determined safely, ask for clarification using `vscode_askQuestions` instead of making assumptions.
- If `copilot-office/<mission-name>/copilot-active-plan.md` is missing, stale, or no longer matches the active goal, run a custom agent named exactly 'Plan' as subagent to produce updated active-plan content using the `## Active Plan Rules` section below, then write that content to `copilot-office/<mission-name>/copilot-active-plan.md` in the main agent.
- In that active-plan generation step, break the active goal into several high-level slices or expected stages when useful. Keep those slices strategic and lightweight rather than turning the active plan into an execution checklist.
- Do NOT refresh the active plan on every loop by default. Refresh it only when the active goal, stage slicing, blocker model, key decisions, or overall strategy changed materially, or when the active plan is missing, stale, or inconsistent.
- Do not move to step 1 until the active goal is resolved safely and the active plan is usable, or a genuine blocker is reported.

## 1. Plan The Next Stage
1. Run a custom agent named exactly 'Plan' as subagent to update or generate the next stage plan:
  - If the `Plan` subagent is unavailable, stop and report that blocker instead of silently substituting another agent.
  1. Confirm the resolved mission folder and ignore unrelated mission folders unless the user explicitly asks for cross-mission work.
  2. Read the mission context from:
    - `copilot-office/<mission-name>/copilot-project-plan.md`
    - `copilot-office/<mission-name>/copilot-active-plan.md`
    - `copilot-office/<mission-name>/copilot-stage-plan.md`
    - `copilot-office/<mission-name>/copilot-desk/`
    - `copilot-office/codebase/CODEBASE.md` when shared architecture matters
  3. Use the slices or expected stages already recorded in `copilot-office/<mission-name>/copilot-active-plan.md` when they are still valid. Re-slice only if the active plan is missing that structure or it is no longer valid.
  4. Produce the next stage-plan content for `copilot-office/<mission-name>/copilot-stage-plan.md` and return it to the main agent.
  5. Keep `copilot-office/<mission-name>/copilot-stage-plan.md` detailed, concrete, and execution-oriented:
    - Capture the current stage objective, ordered tasks, dependencies, verification steps, and immediate blockers.
    - Replace or rewrite this file when the current stage is completed or re-scoped rather than appending a long history.

2. Update `copilot-office/<mission-name>/copilot-stage-plan.md` in the main agent with the stage-plan content returned by the `Plan` subagent.


  


## 2. Execute The Current Stage
Run a subagent to execute the plan defined in `copilot-office/<mission-name>/copilot-stage-plan.md`:
- When delegating to subagents, include these:
  1. the resolved mission name
  2. `../instructions/project_context.instructions.md` (to understand project context)
  3. the stage plan file (direct order for subagent)
- Tell them to stay within that mission folder and its related shared docs.
- Use a dedicated subagent for each distinct multi-step task (e.g., implementation, research, testing, debugging).
- Only pause execution if you encounter a critical blocker, require clarification, or hit a defined milestone requiring manual verification.
- Do not silently switch to another mission folder mid-run.
- The execution subagent owns substantive implementation, research, testing, and debugging work for this stage.
- The main agent may use direct tool calls only for minimal coordination work, such as checking the updated stage files, reading a concise result, or validating whether another handoff is needed.
- Do not use the main agent for substantive stage execution when that work belongs to the execution subagent.
- Do not pass unrelated mission files, the full mission history, the entire `copilot-desk/`, or long raw logs unless the current stage truly depends on them.
- After completing the current stage, reassess the active goal.
- Update `copilot-office/<mission-name>/copilot-active-plan.md` only when it truly needs modification.
- Run the custom agent named exactly `Plan` as subagent to produce updated active-plan content only when the active goal changed, the expected stage slices changed, the blocker model changed materially, key decisions changed, or the current strategy is no longer valid. Follow the `## Active Plan Rules` section below in that active-plan generation step.
- If the active plan only needs a small factual synchronization that does not change goal structure, stage slicing, blockers, key decisions, or strategy, the main agent may update `copilot-office/<mission-name>/copilot-active-plan.md` directly without invoking `Plan`.
- If the active goal is not yet resolved and no genuine blocker exists, return to step 1 with updated state and plan the next stage. Before returning to step 1 for another loop, read `./work_for_goal.prompt.md` again to refresh this workflow and reduce drift across repeated iterations.
- Keep iterating stages until the active goal is resolved or blocked.

## 3. Reporting
Once you reach a natural stopping point, Run a subagent to provide a comprehensive report detailing the resolved mission folder, whether it came from current session history or the inspection step, the mission goal, the active goal, completed stages, remaining backlog, any plan cleanup performed, and any uncovered issues or blockers.

## Active Plan Rules
Keep `copilot-office/<mission-name>/copilot-active-plan.md` actionable, current, and scoped to the active goal:
    - Preserve enough context for a future agent to resume work safely.
    - Keep the active goal, active backlog, blockers, key unresolved decisions, and the minimum recent progress needed to continue.
    - Break the active goal into several high-level slices or expected stages when useful, but you don't need to write details of the slices. Each stage plans will be written separately in `copilot-office/<mission-name>/copilot-stage-plan.md` when that stage is active.
    - Remove or compress irrelevant, outdated, duplicated, or archival content when it is causing confusion.
    - Move durable design notes, reference material, and other non-operational content into `copilot-office/<mission-name>/copilot-desk/` when appropriate, and leave a short pointer in `copilot-office/<mission-name>/copilot-active-plan.md`.