---
name: load_history
description: Load active mission history, including plans, chat logs, and architectural notes.
---

You MUST perform the following context-gathering steps.

## 1. Context Acquisition
1. Resolve `copilot-office/<mission-name>/` using the rules in `../instructions/project_context.instructions.md`. Prefer the mission already established in the current session history. If none exists, inspect the current request and closely related mission references. Ask only if the mission still cannot be resolved safely.
2. Read and retrieve context from:
- `copilot-office/<mission-name>/copilot-project-plan.md`: Digest mission goals.
- `copilot-office/<mission-name>/copilot-current-plan.md`: Digest current actionable tasks.
- `copilot-office/<mission-name>/copilot-chat-logs/`: Read logs in chronological order and summarize them internally to prevent context bloat.
- `copilot-office/<mission-name>/copilot-desk/`: Read files to understand past decisions and reference materials. Summarize heavily.
- `copilot-office/codebase/CODEBASE.md`: Understand shared workspace architecture and conventions that may affect the active mission.
3. Ignore unrelated mission folders unless the user explicitly asks for cross-mission comparison or migration work.

## 2. Status Report
After gathering context, you MUST present a concise and informative status report to the user summarizing:
- The resolved mission folder and the current state of that mission.
- Whether that mission came from current session history or the inspection step.
- The immediate next steps.
- Any conflicting information between mission-specific docs and shared docs.
- Any missing context that requires the user's clarification before proceeding.