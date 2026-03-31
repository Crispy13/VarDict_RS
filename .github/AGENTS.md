# Agent Team — VarDictJava → Rust Parity

This workspace uses a focused parity agent team for VarDictJava-to-Rust translation work. The team is intentionally split by responsibility so Java analysis, Rust implementation, parity validation, and code review stay specialized, while a single orchestrator coordinates delegation and progress.

| Agent Name | File | Role | User-Invocable | Tool Categories |
|---|---|---|---|---|
| Parity Orchestrator | parity-orchestrator.agent.md | Lead coordinator. Delegates to all other agents. Does NOT write code or run commands. | yes | read, search, agent delegation, todo |
| Java Analyst | java-analyst.agent.md | Read-only analysis of VarDictJava source code. | no | read, search, web |
| Rust Implementer | rust-implementer.agent.md | Writes Rust code for parity. | no | read, search, edit, execute, web |
| Parity Tester | parity-tester.agent.md | Validates byte-identical output. | no | read, search, edit, execute |
| Code Reviewer | code-reviewer.agent.md | Reviews ported code quality. | no | read, search, execute |

## Workflow

The Parity Orchestrator delegates work to four specialists:

`Parity Orchestrator -> {Java Analyst, Rust Implementer, Parity Tester, Code Reviewer}`

In practice, the orchestrator scopes the task, requests Java-side analysis, hands implementation work to the Rust specialist, asks the parity tester to validate byte-identical behavior, and uses the reviewer as the final quality gate before closing the loop.

## Other Modes

`🌯GeneralDirector` mode exists for general-purpose tasks in this workspace, but it is not part of the VarDictJava parity team documented here. Its internal structure and subagent setup are intentionally not documented in this file.