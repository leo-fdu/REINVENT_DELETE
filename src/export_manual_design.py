"""Export saved manual designs as per-target REINVENT input files."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
MODES = ("libinvent", "linkinvent")


def export(designs: Path, output: Path) -> dict:
    raw = designs.read_bytes()
    snapshot = json.loads(raw)
    if snapshot.get("schema_version") != 2:
        raise ValueError("Expected manual design schema_version 2")
    if output.exists():
        raise ValueError(f"Output already exists; choose a new directory: {output}")

    files = {}
    tasks = []
    counts = Counter()
    for target, info in sorted(snapshot["targets"].items()):
        if not re.fullmatch(r"[a-zA-Z0-9_-]+", target):
            raise ValueError(f"Invalid target name: {target}")
        if set(info["tasks"]) != set(MODES):
            raise ValueError(f"Expected both LibINVENT and LinkINVENT tasks: {target}")
        source = ROOT / info["source_mol2"]
        if hashlib.sha256(source.read_bytes()).hexdigest() != info["source_sha256"]:
            raise ValueError(f"Original ligand has changed: {target}")
        for mode in MODES:
            task = info["tasks"][mode]
            status = task["status"]
            if status not in ("designed", "skipped"):
                raise ValueError(f"Task is not finalized: {target}/{mode} ({status})")
            counts[status] += 1
            record = {
                "target": target,
                "mode": mode,
                "status": status,
                "smiles_file": None,
                "source_mol2": info["source_mol2"],
                "source_sha256": info["source_sha256"],
                "note": task.get("note", ""),
            }
            if status == "designed":
                smiles = task["input_smiles"]
                if (
                    not isinstance(smiles, str)
                    or not smiles
                    or any(c.isspace() for c in smiles)
                ):
                    raise ValueError(f"Expected one nonempty SMILES line: {target}/{mode}")
                parts = smiles.split("|")
                expected_parts = 1 if mode == "libinvent" else 2
                if len(parts) != expected_parts or any("*" not in part for part in parts):
                    raise ValueError(f"Unexpected input fragment format: {target}/{mode}")
                if task["validation"].get("connectivity_match") is not True:
                    raise ValueError(f"Design connectivity did not match: {target}/{mode}")
                relative = f"{target}/{mode}.smi"
                files[relative] = smiles + "\n"
                record.update(
                    smiles_file=relative,
                    input_smiles=smiles,
                    validation=task["validation"],
                )
            tasks.append(record)

    if not files:
        raise ValueError("No designed tasks to export")
    actual_counts = {
        status: counts[status]
        for status in ("designed", "skipped", "incomplete", "not_started")
    }
    if actual_counts != snapshot["status_counts"]:
        raise ValueError("Snapshot status_counts do not match its tasks")
    resolved_source = designs.resolve()
    source_label = (
        str(resolved_source.relative_to(ROOT))
        if resolved_source.is_relative_to(ROOT)
        else str(resolved_source)
    )
    manifest = {
        "schema_version": 1,
        "source_designs": source_label,
        "source_designs_sha256": hashlib.sha256(raw).hexdigest(),
        "source_created_at": snapshot["created_at"],
        "rdkit_version": snapshot["rdkit_version"],
        "model_ids": snapshot["model_ids"],
        "status_counts": actual_counts,
        "tasks": tasks,
    }
    output.mkdir(parents=True)
    for relative, text in files.items():
        path = output / relative
        path.parent.mkdir(exist_ok=True)
        path.write_text(text, encoding="utf-8")
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--designs", type=Path, required=True, help="Saved manual designs.json")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    try:
        manifest = export(args.designs, args.output)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"Export failed: {error}\n")
    counts = manifest["status_counts"]
    print(
        f"Exported {counts['designed']} input files; "
        f"skipped {counts['skipped']} tasks: {args.output}"
    )


if __name__ == "__main__":
    main()
