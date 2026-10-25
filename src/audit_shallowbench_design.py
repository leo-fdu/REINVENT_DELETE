"""Enumerate allowed cuts to audit compatibility, without choosing user tasks."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import json
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "real-world_dataset/task_design_manual"))
import design


def main():
    design.configure_dataset(ROOT / "shallowbench-7targets/task_design/dataset.json")
    priors = design.models()
    RDLogger.DisableLog("rdApp.warning")
    report = {"scope": "All non-ring single cuts in the editor; compatibility is not biological suitability",
              "model_ids": {mode: prior["metadata"]["model_id"] for mode, prior in priors.items()}, "ligands": {}}
    for key, info in design.atlas.entries.items():
        mol = design.ligand(key)
        bonds = sorted(design.eligible_bonds(mol))
        item = {"heavy_atoms": mol.GetNumHeavyAtoms(), "mw": Descriptors.MolWt(mol),
                "formal_charge": Chem.GetFormalCharge(mol), "smiles": Chem.MolToSmiles(mol),
                "editable": info["editable"], "eligible_bonds": bonds,
                "center_distance_angstrom": info["center_distance_angstrom"], "modes": {}}
        report["ligands"][key] = item
        if not info["editable"]:
            continue
        for mode in design.MODES:
            counts, failures, examples = Counter(), Counter(), []
            for cuts in combinations(bonds, 1 if mode=="libinvent" else 2):
                parts = design.fragments_for(key, list(cuts))
                retained_options = [[p] for p in parts] if mode=="libinvent" else [
                    [p for p in parts if len(design.atlas.dummies(p))==1]]
                for kept in retained_options:
                    spec = {"cuts": [list(c) for c in cuts], "retained": [sorted(design.atlas.ids(p)) for p in kept]}
                    counts["attempts"] += 1
                    try:
                        result = design.build_task(key, mode, spec, priors[mode])
                    except (ValueError, AssertionError) as error:
                        counts["rejected"] += 1
                        failures[str(error)] += 1
                        continue
                    counts["accepted"] += 1
                    stereo = result["validation"]["stereochemistry_match"]
                    counts["stereo_exact" if stereo else "stereo_mismatch"] += 1
                    sizes = [len(p["atom_ids"]) for p in result["retained"]]
                    if min(sizes) < 3:
                        counts["retained_fragment_under_3_heavy_atoms"] += 1
                    examples.append({"cuts": spec["cuts"], "retained": spec["retained"],
                                     "retained_heavy_atoms": sizes,
                                     "generated_heavy_atoms": len(result["generated_reference"]["atom_ids"]),
                                     "input_tokens": result["input_tokens"],
                                     "stereochemistry_match": stereo,
                                     "input_smiles": result["input_smiles"]})
            item["modes"][mode] = {"counts": dict(counts), "failures": dict(failures), "compatible_cuts": examples}
        print(key, {mode: result["counts"] for mode,result in item["modes"].items()}, flush=True)
    out = ROOT / "shallowbench-7targets/task_design/compatibility_audit.json"
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")


if __name__ == "__main__":
    main()
