"""Load one local PLANET model and precompute one protein pocket per process."""

import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import sys

from .chemistry import checked_pocket, prepare_molecule


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PlanetPredictor:
    def __init__(self, cfg):
        import numpy as np
        import torch

        self.torch = torch
        self.cfg = cfg
        device = torch.device(cfg["device"])
        if device.type not in ("cpu", "cuda"):
            raise ValueError("Only cpu and cuda devices are supported")
        if device.type == "cuda":
            if not torch.cuda.is_available():
                raise ValueError(f"Requested CUDA device unavailable: {device}")
            index = device.index if device.index is not None else torch.cuda.current_device()
            if index >= torch.cuda.device_count():
                raise ValueError(f"Requested CUDA device unavailable: {device}")
            device = torch.device(f"cuda:{index}")
        self.device = str(device)
        random.seed(cfg["seed"])
        np.random.seed(cfg["seed"])
        torch.manual_seed(cfg["seed"])
        if device.type == "cuda":
            torch.cuda.manual_seed_all(cfg["seed"])
        root = Path(cfg["planet_root"])
        sys.path.insert(0, str(root))
        self.upstream = importlib.import_module("chemutils")
        model_module = importlib.import_module("PLANET_model")
        for name in ("chemutils", "PLANET_model", "layers", "nnutils"):
            if Path(sys.modules[name].__file__).resolve().parent != root:
                raise RuntimeError(f"Module {name} was imported from outside planet_root")
        self.model = model_module.PLANET(300, 8, 300, 300, 3, 10, 1, device=device)
        self.model.load_state_dict(torch.load(cfg["checkpoint"], map_location=device, weights_only=True))
        self.model.to(device).eval().requires_grad_(False)
        self.pocket, self.pocket_info = checked_pocket(self.upstream, cfg)
        with torch.inference_mode():
            self.res_features = self.model.cal_res_features_helper(
                self.pocket.res_features, self.pocket.alpha_coordinates)
        if not torch.isfinite(self.res_features).all():
            raise RuntimeError("Non-finite precomputed protein features")
        # Include actual pocket tensors: upstream set/alternate ordering may depend on PYTHONHASHSEED.
        pocket_digest = hashlib.sha256()
        for tensor in (self.pocket.res_features, self.pocket.alpha_coordinates):
            pocket_digest.update(tensor.numpy().tobytes())
        versions = {name: importlib.metadata.version(name)
                    for name in ("torch", "rdkit", "numpy", "scipy", "pandas")}
        self.manifest = {
            "version": 1, "config": cfg, "device": self.device,
            "device_name": torch.cuda.get_device_name(device) if device.type == "cuda" else platform.machine(),
            "versions": {"python": platform.python_version(), **versions},
            "cuda_runtime": torch.version.cuda, "python_hash_seed": os.environ.get("PYTHONHASHSEED"),
            "pocket": self.pocket_info, "pocket_tensor_sha256": pocket_digest.hexdigest(),
            "source_sha256": {p.name: sha256(p) for p in sorted(root.glob("*.py"))},
            "adapter_sha256": {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
            "input_sha256": {key: sha256(cfg[key]) for key in ("checkpoint", "protein_pdb", "ligand_sdf") if key in cfg},
        }
        self.oracle_id = hashlib.sha256(json.dumps(self.manifest, sort_keys=True, allow_nan=False).encode()).hexdigest()
        self.manifest["oracle_id"] = self.oracle_id

    def prepare(self, smiles):
        return prepare_molecule(smiles, self.upstream.MAX_NB)

    def predict(self, molecules):
        with self.torch.inference_mode():
            graph = self.upstream.mol_batch_to_graph(molecules, auto_detect=False)
            features, scope = self.model.cal_res_features(self.res_features, len(molecules))
            return self.model.screening(features, scope, graph).reshape(-1).cpu().tolist()
