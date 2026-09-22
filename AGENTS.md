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
2. Commit the changes.
3. Use a short, descriptive commit message.

Do not leave completed work uncommitted.