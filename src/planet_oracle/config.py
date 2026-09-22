"""Strict JSON configuration; paths are relative to the configuration file."""

import json
import math
from pathlib import Path


def reject_constant(value):
    raise ValueError(f"Non-finite JSON constant: {value}")


def loads(value):
    return json.loads(value, parse_constant=reject_constant)


def read_config(path):
    data = loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("Configuration must be a JSON object")
    return data


def finite_number(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return value


def nonempty_string(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def integer(value, name, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an integer in [{low}, {high}]")
    return value


def server_config(path):
    cfg = read_config(path)
    allowed = {"target_id", "planet_root", "checkpoint", "protein_pdb", "ligand_sdf",
               "center", "device", "batch_size", "port", "seed"}
    if cfg.keys() - allowed:
        raise ValueError(f"Unknown server configuration: {sorted(cfg.keys() - allowed)}")
    nonempty_string(cfg.get("target_id"), "target_id")
    base = Path(path).resolve().parent

    def resolve(value, name, directory=False):
        result = (base / nonempty_string(value, name)).resolve()
        if not (result.is_dir() if directory else result.is_file()):
            raise ValueError(f"Missing {name}: {result}")
        return str(result)

    cfg["planet_root"] = resolve(cfg.get("planet_root"), "planet_root", directory=True)
    cfg["protein_pdb"] = resolve(cfg.get("protein_pdb"), "protein_pdb")
    cfg["checkpoint"] = resolve(
        cfg.get("checkpoint", str(Path(cfg["planet_root"]) / "PLANET.param")), "checkpoint")
    if ("ligand_sdf" in cfg) == ("center" in cfg):
        raise ValueError("Specify exactly one of ligand_sdf or center")
    if "ligand_sdf" in cfg:
        cfg["ligand_sdf"] = resolve(cfg["ligand_sdf"], "ligand_sdf")
    else:
        center = cfg["center"]
        if not isinstance(center, list) or len(center) != 3:
            raise ValueError("center must contain three finite numbers")
        for value in center:
            finite_number(value, "center")
    cfg["device"] = nonempty_string(cfg.get("device", "cpu"), "device")
    cfg["batch_size"] = integer(cfg.get("batch_size", 32), "batch_size", 1, 2**31 - 1)
    cfg["port"] = integer(cfg.get("port", 8765), "port", 1, 65535)
    cfg["seed"] = integer(cfg.get("seed", 42), "seed", 0, 2**32 - 1)
    return cfg
