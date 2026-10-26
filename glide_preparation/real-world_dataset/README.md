# MOE preprocessing inputs for six REINVENT cohorts

These are SMILES inputs, not prepared 3D ligands or Glide docking inputs.
Raw results are unchanged. Each of six run folders contains 30 headerless
`<target>_<mode>.smi` files, one `manifest.csv`, and one `summary.json`.
The six runs cover the same 15 targets and LibINVENT/LinkINVENT tasks; ACES
is absent from these generation results. Do not manufacture an ACES task.

## Selection and chemistry

- `reinvent-20260927_first1000`: first 1,000 source records per task, in CSV
  occurrence order, from the original 10,000-attempt tasks.
- The other five folders use all 1,000 source records per task from their
  independent runs. The corresponding source folder name is unchanged.
- Extract the complete molecule from `smiles`, not `generated_smiles` or
  `input_smiles`. Preserve the original SMILES, stereochemistry and duplicates.
- No ranking, deduplication, state enumeration, neutralization, salt removal,
  3D generation or optimization is performed. No replacement of rejected rows.
- Exclude source-invalid records (`smiles_state == 0`), empty/whitespace-bearing
  SMILES, RDKit parsing failures, dummy atoms, and disconnected structures from
  `.smi`. Record every selected attempt and all applicable rejection reasons in
  the manifest. RDKit sanitization is a basic input check, not MOE acceptance.

Each `.smi` line is `SMILES<TAB>molecule_id` with no header. During MOE import,
interpret column 1 as SMILES and column 2 as the molecule name/ID. Keep this ID
as the title or a property in exported SDF; check that it is not truncated.
If preparation creates multiple states, retain the original ID as `parent_id`
and assign distinct state IDs. Run names here identify cohorts, not random seeds.

`source_data_row`, `task_record`, and `smi_line` are 1-based; source_data_row
excludes the CSV header. Rejected rows have an empty smi_line. Source metadata
(seed, step, sample_index, smiles_state) is retained. summary.json records the
compressed source SHA-256, RDKit version and per-task counts/rejection reasons.

Each `.smi` also has a same-named `.txt` containing only one SMILES per line,
without a header, ID or annotation column, as requested for MOE import.
The TXT and SMI records have identical order and counts. Keep the manifests
for tracing TXT line numbers back to source records using `smi_line`.

| Run | Selected | Exported | Rejected |
|---|---:|---:|---:|
| Original first 1,000 | 30,000 | 29,380 | 620 |
| high8 | 30,000 | 29,637 | 363 |
| high9 | 30,000 | 29,608 | 392 |
| high9 cap9 size .01 | 30,000 | 29,722 | 278 |
| high9 cap9 size .015 | 30,000 | 29,753 | 247 |
| high9 cap9 size .02 | 30,000 | 29,768 | 232 |
| Total | 180,000 | 177,868 | 2,132 |

## Reproduce

From the repository root, with Python and RDKit installed:

```bash
python glide_preparation/export_smiles.py
```

The exporter validates source budgets and matching task coverage before writing
each cohort. Re-running replaces the named derived files deterministically;
keep manual MOE outputs elsewhere. Generation used RDKit 2026.03.5 for this
export validation. This does not claim MOE or Glide execution or validation.

## Next handoff

Import a small `.smi` file/sample into MOE first. Confirm names, molecule counts,
3D generation, protonation policy, stereochemistry and minimization settings.
Export prepared 3D SDF with explicit hydrogens and retained IDs, and check its
import into Maestro before batch processing. Store rejected preparation records.
Glide docking also requires a separately generated receptor grid and a saved
SP job configuration. Existing Vina PDBQT/grid JSON files are not Glide grids.
