---
description: Default behavior rules and constraints for this workspace.
applyTo: '**'
---

## 1. Environment and Tools
- **Python / Build Environment:** You MUST ALWAYS prioritize using the `rust_build_env` conda environment (`conda activate rust_build_env`). There is no need to source `conda.sh` unless you invoke a new shell. 
- **Conda/venv Fallback:** If the `rust_build_env` Conda environment does not work, DO NOT guess or fail silently. IMMEDIATELY stop, report the exact problem to the user, and explicitly ask whether you should continue using Conda or switch to a Python `venv`.
- **Temporary Files:** You MUST ALWAYS use `./tmp` for creating any temporary or intermediate files. NEVER use the system `/tmp` directory.

## 2. Problem-Solving Methodology
- **When Stuck:** Before explicitly asking the user for help, you MUST attempt to solve the problem systematically yourself. If you remain stuck, explicitly state what information you are missing and what specific contexts you need from the user to proceed.
- **Alternative Approaches:** If you discover a better or more efficient way to solve a problem than what the user requested, ALWAYS suggest it. However, you MUST clearly explain *why* your alternative option is better than the current approach.

## 3. Context Management
- **Terminal Output:** Redirect stdout and stderr to files in `./tmp` when running commands. Use redirection to terminal when you really need it.