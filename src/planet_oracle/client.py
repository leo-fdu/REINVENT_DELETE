"""Standard-library stdin/HTTP/stdout bridge for REINVENT ExternalProcess."""

import argparse
import json
import sys
from urllib.parse import urlsplit
from urllib.request import ProxyHandler, Request, build_opener

from . import PROTOCOL_VERSION
from .config import finite_number, loads, nonempty_string, read_config

# Finite placeholder for REINVENT's numeric endpoint; the configured increasing
# sigmoid maps it to exactly 0.0 while planet_affinity remains null.
INVALID_AFFINITY_SENTINEL = -1_000_000.0


def client_config(path):
    cfg = read_config(path)
    allowed = {"url", "target_id", "oracle_id", "timeout"}
    if cfg.keys() - allowed:
        raise ValueError(f"Unknown client configuration: {sorted(cfg.keys() - allowed)}")
    for key in ("url", "target_id", "oracle_id"):
        nonempty_string(cfg.get(key), key)
    url = urlsplit(cfg["url"])
    if url.scheme not in ("http", "https") or not url.hostname or url.query or url.fragment or url.username or url.password:
        raise ValueError("url must be an HTTP(S) service base URL without credentials/query/fragment")
    cfg["url"] = cfg["url"].rstrip("/")
    cfg["timeout"] = finite_number(cfg.get("timeout", 120), "timeout")
    if cfg["timeout"] <= 0:
        raise ValueError("timeout must be positive")
    return cfg


def validate_response(data, count, cfg):
    if not isinstance(data, dict):
        raise ValueError("Response must be a JSON object")
    if type(data.get("version")) is not int or data["version"] != PROTOCOL_VERSION:
        raise ValueError("Response protocol version mismatch")
    for key in ("target_id", "oracle_id"):
        if data.get(key) != cfg[key]:
            raise ValueError(f"Response {key} mismatch")
    rows = data.get("results")
    if not isinstance(rows, list) or len(rows) != count:
        raise ValueError("Response results length mismatch")
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or set(row) != {"index", "affinity", "status", "error"}:
            raise ValueError("Malformed result record")
        if type(row["index"]) is not int or row["index"] != index:
            raise ValueError("Missing, duplicate or out-of-order result index")
        if row["status"] == "ok":
            finite_number(row["affinity"], "affinity")
            if row["error"] is not None:
                raise ValueError("Successful result contains an error")
        elif row["status"] == "invalid_input":
            if row["affinity"] is not None:
                raise ValueError("Invalid input must not have a predicted affinity")
            nonempty_string(row["error"], "invalid_input error")
        else:
            raise ValueError("Unknown result status")
    return rows


def score(smiles, cfg):
    request = Request(cfg["url"] + "/score", method="POST",
                      headers={"Content-Type": "application/json"},
                      data=json.dumps({"version": PROTOCOL_VERSION, "target_id": cfg["target_id"],
                                       "smiles": smiles}, allow_nan=False).encode("utf-8"))
    # Loopback traffic must not be redirected through environment-configured proxies.
    with build_opener(ProxyHandler({})).open(request, timeout=cfg["timeout"]) as response:
        rows = validate_response(loads(response.read()), len(smiles), cfg)
    return {"version": PROTOCOL_VERSION, "payload": {
        # REINVENT requires one finite numeric value per row before applying its transform.
        # Preserve the true nullable affinity separately for reporting and provenance.
        "planet_affinity_for_scoring": [
            row["affinity"] if row["status"] == "ok" else INVALID_AFFINITY_SENTINEL
            for row in rows
        ],
        "planet_affinity": [row["affinity"] for row in rows],
        "planet_status": [row["status"] for row in rows],
        "planet_error": [row["error"] or "" for row in rows],
        "planet_target_id": [cfg["target_id"]] * len(rows),
        "planet_oracle_id": [cfg["oracle_id"]] * len(rows),
    }}


def stdin_smiles(stream):
    # Iterating over lines retains internal blank lines; one final newline terminates a line.
    return [line.rstrip("\r\n") for line in stream]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    try:
        cfg = client_config(args.config)
        data = score(stdin_smiles(sys.stdin), cfg)
        print(json.dumps(data, allow_nan=False))
    except Exception as exc:
        print(f"PLANET client failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
