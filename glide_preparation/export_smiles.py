"""Export six fixed REINVENT cohorts without changing molecular chemistry."""
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "glide_preparation/real-world_dataset"
RUNS = [
    ("reinvent-20260927", 10000),
    ("reinvent-high8-1000-20261025", 1000),
    ("reinvent-high9-1000-20261025", 1000),
    ("reinvent-high9-cap9-size001-1000-20261025", 1000),
    ("reinvent-high9-cap9-size0015-1000-20261025", 1000),
    ("reinvent-high9-cap9-size002-1000-20261025", 1000),
]


def rejection(row):
    """Report all applicable problems; never repair or replace a source SMILES."""
    reasons = []
    smiles = row["smiles"]
    if row["smiles_state"] == "0":
        reasons.append("source_invalid")
    if not smiles or any(c.isspace() for c in smiles):
        return ";".join(reasons + ["empty_or_whitespace"])
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ";".join(reasons + ["rdkit_parse_failed"])
    if not mol.GetNumAtoms():
        reasons.append("empty_molecule")
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        reasons.append("dummy_atom")
    if len(Chem.GetMolFrags(mol)) != 1:
        reasons.append("multiple_fragments")
    return ";".join(reasons)


def main():
    RDLogger.DisableLog("rdApp.error")
    RDLogger.DisableLog("rdApp.warning")
    expected_tasks = None
    for run, budget in RUNS:
        source = ROOT / "results" / run / "all_generated.csv.gz"
        counts, selected = Counter(), defaultdict(list)
        with gzip.open(source, "rt", newline="") as handle:
            for source_row, row in enumerate(csv.DictReader(handle), 1):
                key = (row["target"], row["mode"])
                counts[key] += 1
                if counts[key] <= 1000:
                    selected[key].append((source_row, row))
        assert len(counts) == 30 and set(counts.values()) == {budget}, counts
        assert {mode for _, mode in counts} == {"libinvent", "linkinvent"}
        if expected_tasks is None:
            expected_tasks = set(counts)
        assert set(counts) == expected_tasks
        name = run + ("_first1000" if budget == 10000 else "")
        folder = DEST / name
        folder.mkdir(parents=True, exist_ok=True)
        fields = ["molecule_id", "target", "method", "mode", "run_id",
                  "source_file", "source_data_row", "task_record", "seed",
                  "step", "sample_index", "smiles_state", "smiles",
                  "export_status", "rejection_reason", "smi_file", "smi_line"]
        summary = []
        with (folder / "manifest.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for (target, mode), rows in sorted(selected.items()):
                filename = f"{target}_{mode}.smi"
                lines, reasons = [], Counter()
                for ordinal, (source_row, row) in enumerate(rows, 1):
                    molecule_id = f"{name}__{target}__{mode}__{ordinal:04d}"
                    reason = rejection(row)
                    if not reason:
                        lines.append(f"{row['smiles']}\t{molecule_id}\n")
                    else:
                        reasons[reason] += 1
                    writer.writerow(dict(
                        molecule_id=molecule_id, target=target, method="reinvent",
                        mode=mode, run_id=name, source_file=str(source.relative_to(ROOT)),
                        source_data_row=source_row, task_record=ordinal,
                        **{key: row[key] for key in ("seed", "step", "sample_index", "smiles_state", "smiles")},
                        export_status="rejected" if reason else "exported",
                        rejection_reason=reason, smi_file=filename,
                        smi_line="" if reason else len(lines)))
                (folder / filename).write_text("".join(lines))
                (folder / filename).with_suffix(".txt").write_text(
                    "".join(line.split("\t", 1)[0] + "\n" for line in lines))
                summary.append(dict(target=target, mode=mode, source_records=counts[target, mode],
                                    selected=len(rows), exported=len(lines),
                                    rejected=len(rows)-len(lines), rejection_reasons=dict(reasons)))
        report = dict(source_file=str(source.relative_to(ROOT)),
                      source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      rdkit_version=rdBase.rdkitVersion, tasks=summary,
                      selected=sum(x["selected"] for x in summary),
                      exported=sum(x["exported"] for x in summary),
                      rejected=sum(x["rejected"] for x in summary))
        (folder / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        print(name, report["selected"], report["exported"], report["rejected"])


if __name__ == "__main__":
    main()
