---
name: work_for_request
description: Interactive plan-work-report cycle for handling isolated user requests or side quests.
---

You MUST perform the following interactive workflow:

1. **Scope the Request:** Decide whether the request is related to an active mission under `copilot-office/<mission-name>/` or unrelated to any active mission. A request may be related to a mission without being part of that mission. In that case, you may refer to mission context, but you are not allowed to modify mission files unless the user explicitly asks.
   - If you cannot determine whether the request is related to an active mission, you MUST ask the user for clarification using `vscode_askQuestions` instead of making assumptions.
   - If it is related to an active mission, resolve the mission folder using `../instructions/project_context.instructions.md`, using current session history first and a narrow inspection step second, and read only the relevant mission files for context.
   - If it has nothing to do with any active mission, do not read or update any mission planning files unless the user explicitly asks.
   - Ask the user with `vscode_askQuestions` only if the mission still cannot be resolved safely.

2. **Interactive Planning:** Use skill 'plan' to generate actionable steps resolving the user's request and present the plan directly within your chat message.
   - For side quests, the plan should be generated from the user's request and the history of the current session rather than from mission planning files.
   - You MUST use the `vscode_askQuestions` tool to seek the user's explicit confirmation.
   - Do NOT place the plan text inside the tool's question property; output the plan plainly in the main chat.
   - If the user rejects the plan, you MUST restart Step 1 incorporating their feedback.

3. **Execution Phase**
You MUST implement the tasks defined in the plan.
- Use a dedicated subagent for each distinct task (e.g., implementation, research, testing, debugging).
- Proceed sequentially through the backlog WITHOUT stopping between tasks.
- Only pause execution if you encounter a critical blocker, require clarification, or hit a defined milestone requiring manual verification.

4. **Status Reporting:** Generate a detailed report identifying what was achieved, what the next steps are, whether mission files were intentionally used or intentionally ignored, and whether the mission came from current session history or the inspection step.