"""Resolve traceable protein/pocket corrections without changing design ligands.

Only standard-library code is used; derived receptors belong to the run directory.
"""

import hashlib
import json
import math
from pathlib import Path
import re


def file_sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def pdb_lines(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise ValueError(f"Multiple PDB models are not supported: {path}")
    if any(line.startswith(("ATOM  ", "HETATM")) and len(line) < 54 for line in lines):
        raise ValueError(f"Truncated PDB atom record: {path}")
    return lines


def select_alternates(lines):
    """One alternate per residue: highest mean occupancy, then lexical label."""
    groups = {}
    for line in lines:
        groups.setdefault((line[21], line[22:27]), []).append(line)
    selected = []
    for contents in groups.values():
        labels = sorted({line[16] for line in contents if line[16] != " "})
        if labels:
            def occupancy(label):
                values = [float(line[54:60]) for line in contents if line[16] == label]
                if not all(math.isfinite(v) and 0 <= v <= 1 for v in values):
                    raise ValueError("Invalid alternate occupancy")
                return math.fsum(values) / len(values)
            chosen = min(labels, key=lambda label: (-occupancy(label), label))
            contents = [line for line in contents if line[16] in (" ", chosen)]
        selected.extend(contents)
    return selected


def receptor_chain(path, chain):
    lines = pdb_lines(path)
    # Keep modified polymer residues as HETATM, with their original chemistry.
    # PLANET reads ATOM records only; no modified residue is relabeled as standard.
    modified = {(line[16], line[18:23]) for line in lines if line.startswith("MODRES")}
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM")) and line[21] == chain and (
        line.startswith("ATOM  ") or
        (line.startswith("HETATM") and (line[21], line[22:27]) in modified)
    )]
    if not any(line.startswith("ATOM  ") for line in atoms):
        raise ValueError(f"No protein ATOM records for chain {chain}: {path}")
    atoms = select_alternates(atoms)
    for line in atoms:
        coordinates = [float(line[start:start + 8]) for start in (30, 38, 46)]
        if not all(map(math.isfinite, coordinates)):
            raise ValueError(f"Non-finite receptor coordinates: {path}")
    header = ["REMARK   Derived PLANET receptor; coordinates unchanged.",
              f"REMARK   Source chain {chain}; one occupancy-selected alternate per residue."]
    header += [line for line in lines if line.startswith("CRYST1") or
               (line.startswith("MODRES") and line[16] == chain)]
    return "\n".join(header + atoms + ["TER", "END"]) + "\n"


def reference_coords(path, chain, resname, residue_id):
    """Heavy-atom coordinates of the pocket-reference residue, one per atom."""
    atoms = [line for line in pdb_lines(path) if line.startswith(("ATOM  ", "HETATM"))
             and line[21] == chain and line[17:20].strip() == resname
             and line[22:27].strip() == residue_id]
    atoms = select_alternates(atoms)
    coordinates, names = [], set()
    for line in atoms:
        element = line[76:78].strip()
        if not element:
            raise ValueError(f"Reference atom has no element: {path}")
        if element in ("H", "D", "T"):
            continue
        name = line[12:16].strip()
        xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
        if name in names or not all(map(math.isfinite, xyz)):
            raise ValueError(f"Duplicate/non-finite reference atoms: {path}")
        names.add(name)
        coordinates.append(xyz)
    if not coordinates:
        raise ValueError(f"Missing reference residue {chain}/{resname}/{residue_id}: {path}")
    return coordinates


def reference_center(path, chain, resname, residue_id):
    coordinates = reference_coords(path, chain, resname, residue_id)
    return [math.fsum(x[axis] for x in coordinates) / len(coordinates) for axis in range(3)]


def corrected_inputs(project, config_path=None):
    """Validate source hashes and return resolved corrections plus frozen files."""
    project = Path(project).resolve()
    default = project / "configs/manual_planet_rl/target_inputs.json"
    if config_path is None and not default.exists():
        return {}, []  # Preserve the original interface for other datasets/fixtures.
    path = Path(config_path).resolve() if config_path is not None else default
    cfg = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(cfg, dict) or set(cfg) != {"schema_version", "targets"}
            or type(cfg["schema_version"]) is not int or cfg["schema_version"] != 1
            or not isinstance(cfg["targets"], dict)):
        raise ValueError("Invalid target inputs schema")
    frozen = [path]
    resolved = {}

    def source(spec, reference=False):
        keys = {"source_pdb", "source_sha256", "chain", "origin"}
        if reference:
            keys |= {"resname", "residue_id"}
        if not isinstance(spec, dict) or set(spec) != keys:
            raise ValueError("Invalid receptor/reference fields")
        if not isinstance(spec["chain"], str) or not re.fullmatch(r"[A-Za-z0-9]", spec["chain"]):
            raise ValueError("A single explicit author chain is required")
        if not isinstance(spec["origin"], str) or not spec["origin"].strip():
            raise ValueError("Source origin is required")
        if not isinstance(spec["source_pdb"], str) or not spec["source_pdb"].strip():
            raise ValueError("Source PDB path is required")
        relative = Path(spec["source_pdb"])
        file = (project / relative).resolve()
        if relative.is_absolute() or not file.is_relative_to(project):
            raise ValueError("Corrected input sources must be project-relative")
        digest = spec["source_sha256"]
        if (not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest)
                or file_sha256(file) != digest):
            raise ValueError(f"Corrected input source changed: {file}")
        if reference and (not isinstance(spec["resname"], str)
                          or not re.fullmatch(r"[A-Za-z0-9]{1,3}", spec["resname"])
                          or not isinstance(spec["residue_id"], str)
                          or not re.fullmatch(r"-?[0-9]+[A-Za-z]?", spec["residue_id"])):
            raise ValueError("Invalid reference residue identity")
        frozen.append(file)
        return {**spec, "source_pdb": str(file)}

    for target, correction in cfg["targets"].items():
        if (not re.fullmatch(r"[a-zA-Z0-9_-]+", target) or not isinstance(correction, dict)
                or set(correction) != {"receptor", "pocket_reference", "reason"}):
            raise ValueError("Invalid target correction")
        if not isinstance(correction["reason"], str) or not correction["reason"].strip():
            raise ValueError("Correction reason is required")
        receptor = source(correction["receptor"])
        reference = source(correction["pocket_reference"], reference=True)
        resolved[target] = {
            "receptor": receptor, "pocket_reference": reference, "reason": correction["reason"],
            "center": reference_center(reference["source_pdb"], reference["chain"],
                                       reference["resname"], reference["residue_id"]),
            "receptor_text": receptor_chain(receptor["source_pdb"], receptor["chain"]),
        }
    return resolved, list(dict.fromkeys(frozen))
