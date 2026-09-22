"""Linux real-model acceptance test against upstream VS_SMI_Dataset inference."""

import argparse
import json
from pathlib import Path
import threading

from .client import score
from .config import server_config
from .server import Oracle, make_server


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--planet-root", type=Path, required=True)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory, must not exist")
    args = parser.parse_args()
    from .backend import PlanetPredictor
    from rdkit import Chem
    import numpy as np
    import torch

    args.output_dir.mkdir(parents=True, exist_ok=False)
    root = args.planet_root.resolve()
    cfg_file = args.output_dir / "server.json"
    cfg_file.write_text(json.dumps({
        "target_id": "planet_demo", "planet_root": str(root),
        "protein_pdb": str(root / "demo/adrb2.pdb"),
        "ligand_sdf": str(root / "demo/adrb2_ligand.sdf"),
        "device": args.device, "batch_size": 2, "seed": 42,
    }, indent=2) + "\n")
    cfg = server_config(cfg_file)
    predictor = PlanetPredictor(cfg)
    (args.output_dir / "manifest.json").write_text(json.dumps(predictor.manifest, indent=2) + "\n")
    smiles = []
    for mol in Chem.SDMolSupplier(str(root / "demo/mols.sdf")):
        if mol is not None:
            smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        if len(smiles) == 5:
            break
    if len(smiles) != 5:
        raise RuntimeError("Need five valid molecules in PLANET/demo/mols.sdf")
    smiles.append(smiles[0])
    smi_file = args.output_dir / "demo.smi"
    smi_file.write_text("".join(f"{smi} demo_{i}\n" for i, smi in enumerate(smiles)))

    # Independent official estimator and official SMILES dataset, on the requested device.
    from PLANET_run import PlanetEstimator, VS_SMI_Dataset

    official = PlanetEstimator(torch.device(args.device))
    official.model.to(torch.device(args.device))
    if next(official.model.parameters()).device != torch.device(predictor.device):
        raise AssertionError("Official model device mismatch")
    official.set_pocket_from_ligand(cfg["protein_pdb"], cfg["ligand_sdf"])
    reference = []
    dataset = VS_SMI_Dataset(str(smi_file), batch_size=2)
    with torch.inference_mode():
        for idx in range(len(dataset)):
            graph, batch_smiles, _ = dataset[idx]
            features, scope = official.model.cal_res_features(official.res_features, len(batch_smiles))
            reference.extend(official.model.screening(features, scope, graph).reshape(-1).cpu().tolist())
    if not np.isfinite(reference).all():
        raise AssertionError("Official model produced non-finite values")

    oracle = Oracle(predictor, cfg["target_id"], cfg["batch_size"])
    identities = (id(predictor.model), id(predictor.pocket), id(predictor.res_features))
    with make_server(oracle, port=0) as server:
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()
        client = {"url": f"http://127.0.0.1:{server.server_port}", "target_id": cfg["target_id"],
                  "oracle_id": predictor.oracle_id, "timeout": 120, "low": 4, "high": 10, "k": 0.5}
        try:
            first = score(smiles, client)["payload"]
            second = score(smiles, client)["payload"]
            singleton = [score([smi], client)["payload"]["planet_affinity"][0] for smi in smiles]
            mixed = score([smiles[0], "invalid", "", "C.C", "C*", smiles[0]], client)["payload"]
            assert score([], client)["payload"]["planet_reward"] == []
            assert mixed["planet_status"] == ["ok"] + ["invalid_input"] * 4 + ["ok"]
            assert mixed["planet_reward"][1:5] == [0.0] * 4
            assert mixed["planet_affinity"][1:5] == [None] * 4
        finally:
            server.shutdown()
            worker.join()
    for values in (first["planet_affinity"], second["planet_affinity"], singleton):
        np.testing.assert_allclose(values, reference, atol=1e-4, rtol=1e-4)
    assert identities == (id(predictor.model), id(predictor.pocket), id(predictor.res_features))
    # Closed listener must fail, never synthesize a reward.
    try:
        score(smiles, client)
    except OSError:
        pass
    else:
        raise AssertionError("Scoring succeeded after service shutdown")
    report = {"passed": True, "device": predictor.device, "oracle_id": predictor.oracle_id,
              "atol": 1e-4, "rtol": 1e-4, "smiles": smiles,
              "official": reference, "adapter": first["planet_affinity"], "singleton": singleton}
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, allow_nan=False))


if __name__ == "__main__":
    main()
