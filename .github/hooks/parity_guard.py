#!/usr/bin/env python3

import json
import sys
from pathlib import Path


STATE_PATH = Path("tmp/copilot-parity-guard-state.json")

VALIDATION_MARKERS = (
    "cargo test",
    "cargo bench",
    "integration_test",
    "_fixture_test",
    "parity",
    "validate",
    "validation",
    "benchmark",
    "bench",
    "v a r d i c t".replace(" ", ""),
)

EDIT_TOOL_MARKERS = (
    "apply_patch",
    "create_file",
    "edit_notebook_file",
    "vscode_renamesymbol",
)

SENSITIVE_PATH_MARKERS = (
    "src/bin/vardict.rs",
    "src/mods/",
    "src/data/",
    "src/scopedata/",
    "src/variants/",
)


def load_state() -> dict:
    if not STATE_PATH.exists():
        return {
            "prompt_has_validation": False,
            "validation_seen": False,
            "prompted_for_edit": False,
        }

    try:
        return json.loads(STATE_PATH.read_text())
    except (json.JSONDecodeError, OSError):
        return {
            "prompt_has_validation": False,
            "validation_seen": False,
            "prompted_for_edit": False,
        }


def save_state(state: dict) -> None:
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps(state))


def contains_any(haystack: str, markers: tuple[str, ...]) -> bool:
    return any(marker in haystack for marker in markers)


def emit_continue(system_message: str | None = None) -> None:
    payload = {"continue": True}
    if system_message:
        payload["systemMessage"] = system_message
    print(json.dumps(payload))


def emit_ask(reason: str, system_message: str) -> None:
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "ask",
                    "permissionDecisionReason": reason,
                },
                "systemMessage": system_message,
            }
        )
    )


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    raw = sys.stdin.read()
    raw_lower = raw.lower()

    state = load_state()
    validation_detected = contains_any(raw_lower, VALIDATION_MARKERS)

    if mode == "user-prompt":
        state = {
            "prompt_has_validation": validation_detected,
            "validation_seen": validation_detected,
            "prompted_for_edit": False,
        }
        save_state(state)
        emit_continue()
        return 0

    if mode == "post-tool":
        if validation_detected:
            state["validation_seen"] = True
            save_state(state)
        emit_continue()
        return 0

    if mode != "pre-tool":
        emit_continue()
        return 0

    is_edit_tool = contains_any(raw_lower, EDIT_TOOL_MARKERS)
    touches_sensitive_path = contains_any(raw_lower, SENSITIVE_PATH_MARKERS)

    if not is_edit_tool or not touches_sensitive_path:
        emit_continue()
        return 0

    if state.get("prompt_has_validation") or state.get("validation_seen"):
        emit_continue()
        return 0

    if state.get("prompted_for_edit"):
        emit_continue(
            "Parity-sensitive Rust edits should end with fixture, parity, or benchmark validation."
        )
        return 0

    state["prompted_for_edit"] = True
    save_state(state)
    emit_ask(
        "Parity-sensitive Rust edit without an explicit validation step in the prompt.",
        "This edit touches parity-sensitive Rust files. Name the narrowest follow-up validation, such as a fixture test, manifest-targeted parity test, or benchmark, before continuing.",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())