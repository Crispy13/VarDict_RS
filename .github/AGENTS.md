# Agent Team — VarDictJava → Rust Parity

This workspace uses a focused parity agent team for VarDictJava-to-Rust translation work. The team is intentionally split by responsibility so Java analysis, Rust implementation, parity validation, review, and documentation cache maintenance stay specialized, while a single orchestrator coordinates delegation and progress.

| Agent Name | File | Role | User-Invocable | Tool Categories |
|---|---|---|---|---|
| Parity Orchestrator | parity-orchestrator.agent.md | Lead coordinator. Delegates to all other agents. Does NOT write code or run commands. | yes | read, search, agent delegation, todo |
| Java Analyst | java-analyst.agent.md | Read-only analysis of VarDictJava source code. | no | read, search, web |
| Rust Implementer | rust-implementer.agent.md | Writes Rust code for parity. | no | read, search, edit, execute, web |
| Parity Tester | parity-tester.agent.md | Validates byte-identical output. | no | read, search, edit, execute |
| Code Reviewer | code-reviewer.agent.md | Reviews ported code quality. | no | read, search, execute |
| Codebase Librarian | codebase-librarian.agent.md | Maintains the Java and Rust codebase documentation cache from agent reports. | no | read, search, edit, web |

## Workflow

The Parity Orchestrator delegates work to five specialists:

`Parity Orchestrator -> {Java Analyst, Rust Implementer, Parity Tester, Code Reviewer, Codebase Librarian}`

In practice, the orchestrator scopes the task, requests Java-side analysis, hands implementation work to the Rust specialist, sends validated reports to the codebase librarian, asks the parity tester to confirm byte-identical behavior, and uses the reviewer as the final quality and performance gate before closing the loop.

## Skill-Agent Relationships

| Skill | Primary Agent | When It Applies | Outcome |
|---|---|---|---|
| `parity-check` | Parity Orchestrator | New module ports and end-to-end parity fixes | Canonical new-module and bug-fix workflow |
| `codebase-doc-manage` | Codebase Librarian | Module orientation, cache updates, audits | Cached module docs stay current |
| `shard-diagnosis` | Rust Implementer | Single-shard failures and column-level mismatches | Precise shard diagnosis before code changes |
| `parity-test-run` | Parity Tester | Running targeted parity cases or shard subsets | Reproducible parity validation |
| `parity-workflow` | Parity Tester | Harness operation, cache usage, shard workflow | Consistent parity infrastructure handling |
| `tiered-config-test` | Parity Orchestrator | Expanding config coverage in stages | Ordered config promotion without sweep sprawl |
| `parity-fix-review` | Rust Implementer | Lightweight gate for shard-scoped fixes | Fast pre-report parity review |
| `change-impact-review` | Code Reviewer | Hot-path changes or pre-merge performance checks | Binding performance verdict |
| `mismatch-triage` | Parity Orchestrator | Classifying mismatch severity and owner | Faster routing to the right module |
| `perf-optimization` | Rust Implementer | Throughput regressions or CPU hot paths | Root-caused performance fixes |
| `mem-optimization` | Rust Implementer | RSS regressions or allocation pressure | Memory-focused optimization workflow |
| `pending-modes` | Parity Orchestrator | Somatic or amplicon parity planning before fixtures are complete | Placeholder workflow for unsupported mode parity |

## Request Routing

| Request Type | Start With | Route | Result |
|---|---|---|---|
| Port a new Java module or method | Parity Orchestrator | `parity-check` → Java Analyst → Rust Implementer → Parity Tester → Code Reviewer → Codebase Librarian | End-to-end parity workflow |
| Fix a known parity mismatch | Parity Orchestrator | `parity-check` Bug Fix Track → specialists as needed | Reproduced, fixed, reviewed mismatch |
| Diagnose a shard failure or column diff | Rust Implementer | `shard-diagnosis` → `parity-fix-review` → Parity Tester | Targeted shard fix loop |
| Run or expand parity tests | Parity Tester | `parity-test-run` and, if needed, `tiered-config-test` or `parity-workflow` | Verified parity coverage |
| Update or audit the codebase cache | Codebase Librarian | `codebase-doc-manage` | Fresh module docs and audit results |
| Evaluate hot-path performance risk | Code Reviewer | `change-impact-review`, escalate to `perf-optimization` or `mem-optimization` if needed | Performance verdict with escalation path |
| Plan somatic or amplicon parity work | Parity Orchestrator | `pending-modes` | Fixture and workflow staging for pending modes |

## Other Modes

`🌯GeneralDirector` mode exists for general-purpose tasks in this workspace, but it is not part of the VarDictJava parity team documented here. Its internal structure and subagent setup are intentionally not documented in this file.