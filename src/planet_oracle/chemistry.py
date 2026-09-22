"""Input checks around (not replacements for) the upstream graph builders."""

from rdkit import Chem


def prepare_molecule(smiles, max_neighbors):
    if not smiles.strip():
        return None, "empty_molecule"
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return None, "invalid_smiles"
    if mol.GetNumAtoms() == 0:
        return None, "empty_molecule"
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        return None, "dummy_atom"
    if len(Chem.GetMolFrags(mol)) != 1:
        return None, "disconnected_structure"
    # Match VS_SMI_Dataset: sanitize and AddHs, without neutralization or 3D embedding.
    mol = Chem.AddHs(mol)
    if any(atom.GetDegree() > max_neighbors for atom in mol.GetAtoms()):
        return None, "graph_capacity_exceeded"
    return mol, None


def checked_pocket(upstream, cfg):
    """Audit every selected residue before ProteinPocket can silently skip it."""
    import numpy as np

    if "ligand_sdf" in cfg:
        supplier = Chem.SDMolSupplier(cfg["ligand_sdf"], sanitize=False)
        ligand = supplier[0] if len(supplier) else None
        if ligand is None or not ligand.GetNumAtoms() or not ligand.GetNumConformers():
            raise ValueError("Reference ligand must have atoms and coordinates")
        coordinates = np.asarray(ligand.GetConformer().GetPositions())
        if not np.isfinite(coordinates).all():
            raise ValueError("Non-finite reference ligand coordinates")
        center = upstream.Mol(ligand).compute_centeroid()
        kwargs = {"ligand_sdf": cfg["ligand_sdf"]}
    else:
        center = np.asarray(cfg["center"], dtype=np.float32)
        kwargs = dict(zip(("centeriod_x", "centeriod_y", "centeriod_z"), cfg["center"]))
    if not np.isfinite(center).all():
        raise ValueError("Non-finite pocket center")

    with open(cfg["protein_pdb"], encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.startswith("ATOM")]
    groups = {}
    selected = set()
    for line in lines:
        key = (line[22:27].strip(), line[21])
        xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
        if not np.isfinite(xyz).all():
            raise ValueError(f"Non-finite PDB coordinates in residue {key}")
        groups.setdefault(key, []).append(line)
        if upstream.near_pocket(line, center):
            selected.add(key)
    if not selected:
        raise ValueError("Empty pocket (12 Angstrom selection)")
    errors = []
    for key in sorted(selected):
        try:
            residue = upstream.Residue(groups[key])
            alpha = residue.get_alpha_position()
            if alpha.shape != (1, 3) or not np.isfinite(alpha).all():
                raise ValueError("expected exactly one finite C-alpha after alternate selection")
            if not np.isfinite(residue.coordinates).all() or not np.isfinite(residue.get_mass_center()).all():
                raise ValueError("non-finite residue coordinates/mass center")
            if not np.isfinite(residue.get_feature()).all():
                raise ValueError("non-finite residue features")
        except Exception as exc:
            errors.append(f"{key}: {type(exc).__name__}: {exc}")
    if errors:
        raise ValueError("Incomplete pocket; upstream would skip or misparse residues: " + "; ".join(errors))
    pocket = upstream.ProteinPocket(cfg["protein_pdb"], **kwargs)
    if pocket.res_count != len(selected):
        raise ValueError("Upstream skipped residues despite pocket audit")
    if tuple(pocket.alpha_coordinates.shape) != (len(selected), 3):
        raise ValueError("Upstream pocket C-alpha count mismatch")
    if not np.isfinite(pocket.alpha_coordinates.numpy()).all() or not np.isfinite(pocket.res_features.numpy()).all():
        raise ValueError("Non-finite upstream pocket tensors")
    return pocket, {"residue_count": pocket.res_count, "radius_angstrom": 12.0,
                    "center": center.tolist(), "residues": [list(key) for key in sorted(selected)]}
