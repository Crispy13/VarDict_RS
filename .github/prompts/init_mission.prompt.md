---
name: init_mission
description: Custom prompt to initialize the project workflow, workspace, and active mission plan.
---

You MUST perform the following initialization steps. Use subagents when useful, but keep the flow centered on one resolved mission folder.

1. **Resolve Mission Folder:** If the user doesn't give mission name, you must ask user it with `vscode_askQuestions` tool. This is init phase, so there must be no session history. Resolve `copilot-office/<mission-name>/` using the rules in `../instructions/project_context.instructions.md`. 

2. **Understand Workspace:** Run subagents to analyze the current workspace and codebase context to understand the requested mission and the surrounding architecture.

3. **Initialize Mission Structure:** Ensure `copilot-office/`, `copilot-office/codebase/`, `copilot-office/<mission-name>/` and `copilot-office/<mission-name>/copilot-desk/` exist.

4. **Locate or Draft Project Plan:** Read `copilot-office/<mission-name>/copilot-project-plan.md` to grasp the mission objective and scope. If the file is missing, draft a minimal project plan from the workspace context and the user request, clearly mark it as provisional, and tell the user that it should be reviewed.

5. **Draft Current Plan:** Read or create `copilot-office/<mission-name>/copilot-current-plan.md` to direct immediate next steps. Keep this plan focused on the active mission only.

6. **Report to User:** Present the resolved mission folder, whether it came from current session history or the inspection step, parsed project plan, immediate next steps, and any missing context concisely. If shared docs appear to conflict with mission-specific docs, call that out explicitly.