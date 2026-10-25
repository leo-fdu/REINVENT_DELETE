"""Capped PLANET sigmoid minus the heavy atoms added to fixed fragments."""

import math


def generated_heavy_atoms(product, fixed_fragments):
    from rdkit import Chem

    def count(smiles):
        mol = Chem.MolFromSmiles(smiles.replace("|", "."))
        if mol is None:
            raise ValueError("Cannot count heavy atoms in an invalid molecule/fragment")
        # Exclude both attachment dummies (atomic number 0) and explicit H.
        return sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())

    added = count(product) - count(fixed_fragments)
    if added < 0:
        raise ValueError("Product has fewer heavy atoms than its retained fragments")
    return added


def reward_terms(affinity, product, fixed_fragments, experiment):
    adjustment = experiment["reward_adjustment"]
    capped = min(float(affinity), adjustment["affinity_cap"])
    if not math.isfinite(capped) or not math.isfinite(float(affinity)):
        raise ValueError("Non-finite PLANET affinity")
    curve = experiment["sigmoid"]
    z = math.log(10) * 10 * curve["k"] / (curve["high"] - curve["low"]) * (
        capped - (curve["high"] + curve["low"]) / 2)
    sigmoid = 1 / (1 + math.exp(-z)) if z >= 0 else math.exp(z) / (1 + math.exp(z))
    added = generated_heavy_atoms(product, fixed_fragments)
    penalty = adjustment["generated_heavy_atom_penalty"] * added
    return {"capped_planet_affinity": capped, "sigmoid_reward": sigmoid,
            "generated_heavy_atoms": added, "size_penalty": penalty,
            "base_reward": sigmoid - penalty}
