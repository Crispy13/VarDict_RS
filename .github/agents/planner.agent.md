---
name: Planner
description: Researches and outlines multi-step plans
argument-hint: Outline the goal or problem to research
target: vscode
user-invocable: false
model: ['Claude Opus 4.6 (copilot)']
tools: ['search', 'read', 'web', 'vscode/memory', 'github/issue_read', 'github.vscode-pull-request-github/issue_fetch', 'github.vscode-pull-request-github/activePullRequest', 'execute/getTerminalOutput', 'execute/testFailure', 'agent', 'vscode/askQuestions']
agents: ['Explore']
handoffs:
  - label: Start Implementation
    agent: agent
    prompt: 'Start implementation'
    send: true
  - label: Open in Editor
    agent: agent
    prompt: '#createFile the plan as is into an untitled file (`untitled:plan-${camelCaseName}.prompt.md` without frontmatter) for further refinement.'
    send: true
    showContinueOn: false
---

Your instruction is here: `$HOME/.vscode-server/data/User/globalStorage/github.copilot-chat/plan-agent/Plan.agent.md`

If you can't resolve the path, then stop and report that you can't find base planning instructions.