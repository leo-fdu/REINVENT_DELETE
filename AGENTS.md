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

## Evaluation Plan

Evaluation lives in `evaluation/` with two independent modules:

```
evaluation/
├── docking/      # AutoDock Vina scoring
└── descriptor/   # molecular descriptors
```

### Data flow

The single, shared contract for both modules is **SMILES strings as input**.

```
model outputs (SMILES / SDF / MOL2)
        │
        ▼
  convert → canonical SMILES table  (target, method, smiles)
        │
        ├──────────────► evaluation/docking/
        └──────────────► evaluation/dcriptor/
```

The docking module is implemented in `evaluation/docking/`
(`prepare_receptors.py` + `run_docking.py`; see its README). Input contract:
one CSV with columns `target,method,smiles` (`method` ∈ `reinvent`/`delete`);
output: `molecule_id,method,target,smiles,vina_score` plus error/summary side
files. Grid boxes and Vina parameters are frozen per target in
`evaluation/docking/configs/*.json`.

DELETE outputs are not SMILES by default; they **must** be converted to SMILES
before entering the evaluation pipeline. The conversion step is a dedicated,
reusable adapter so both DELETE and REINVENT4 produce the same intermediate
format. Raw model outputs are preserved; converted tables are written as new
derived files.

### `evaluation/docking/`

- **Inputs:** generated-molecule SMILES + target 3D structure (receptor).
- **Pipeline:** SMILES → 3D conformer generation (RDKit ETKDG + MMFF/UFF
  minimization) → ligand PDBQT → Vina docking against the prepared receptor
  PDBQT → binding affinity (kcal/mol) per molecule.
- Receptor preparation and the docking grid box are configured per target, one
  config per target shared by both models so the comparison stays fair.
- Outputs: per-molecule scores plus per-target summary statistics (best, mean,
  median, success count).

### `evaluation/descriptor/`

- **Inputs:** the same SMILES table (no 3D structure needed).
- **Pipeline:** SMILES → RDKit descriptors (e.g. MW, LogP, TPSA, HBD, HBA,
  rotatable bonds, ring counts, QED, SA score) and optionally fingerprints.
- Outputs: one descriptor row per molecule, plus per-target summary statistics
  for DELETE vs REINVENT4 comparison.

### Fairness rules

- Both modules consume the identical SMILES table for a given target, so any
  difference in the reported numbers comes from the generation model, not from
  the evaluation path.
- Target structures, grid boxes, and descriptor settings are defined once in
  `configs/` (or `evaluation/` configs) and reused across models.
- Evaluation scripts must be deterministic (fixed random seed for conformer
  generation) and reproducible from raw outputs.

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
