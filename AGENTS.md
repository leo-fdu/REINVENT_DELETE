# AGENTS.md

## Project Goal

This project benchmarks two ligand-generation approaches across 16 protein targets:

- **REINVENT4**
  - LibINVENT
  - LinkINVENT
- **DELETE**
  - Corresponding fragment growing / fragment linking tasks

The goal is to compare the generation performance of the two models under consistent experimental settings across all 16 targets.

## Scope

For REINVENT4, only work on:

- **LibINVENT**
- **LinkINVENT**

Do **not** introduce Mol2Mol or other REINVENT4 generation modes unless explicitly requested.

Experiments should use the same target structures and, as far as possible, equivalent ligand fragments / constraints for both models so that comparisons are fair.

## Development Rules

- Keep scripts, configurations, intermediate data, and results organized and reproducible.
- Do not make unnecessary changes outside the current task.
- Prefer simple and transparent implementations over unnecessary abstractions.
- Preserve raw input data; derived or processed files should be stored separately.

## Git

**Every meaningful change must be committed to Git.**

After completing a coherent unit of work:

1. Check the changes with `git status` / `git diff`.
2. Before committing, audit for files that may be generated, local-only, sensitive, or otherwise inappropriate for version control but are not covered by `.gitignore`.
3. For every possible candidate, report its path, why it may need to be ignored, and the proposed ignore pattern to the user. Ask the user whether the pattern should be added; do not silently add an ignore rule, stage the candidate, or delete it.
4. After the user decides how to handle any candidates, commit the changes.
5. Use a short, descriptive commit message.

The pre-commit audit should include:

```bash
git status --short --untracked-files=all
git status --short --ignored --untracked-files=all
git ls-files --others --exclude-standard
```

Use `git check-ignore -v --no-index -- <path>` to verify whether a specific path is covered by an ignore rule. Review both the root `.gitignore` and the `.gitignore` files belonging to the `Delete` and `REINVENT4` submodules. Do not use broad rules that could hide raw inputs, configurations, reproducible scripts, or intentionally versioned results without the user's decision.

Do not leave completed work uncommitted.
