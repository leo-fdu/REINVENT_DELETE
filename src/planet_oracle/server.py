"""Serial, loopback-only HTTP service. Unexpected inference errors fail the request."""

import argparse
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import logging
from pathlib import Path

from . import PROTOCOL_VERSION
from .config import finite_number, loads, server_config

LOGGER = logging.getLogger(__name__)


class Oracle:
    def __init__(self, predictor, target_id, batch_size=32):
        self.predictor = predictor
        self.target_id = target_id
        self.batch_size = batch_size

    def health(self):
        return {"version": PROTOCOL_VERSION, "ready": True, "target_id": self.target_id,
                "oracle_id": self.predictor.oracle_id, "device": self.predictor.device,
                "pocket": self.predictor.pocket_info}

    def score(self, smiles):
        results, valid = [], []
        for index, smi in enumerate(smiles):
            mol, error = self.predictor.prepare(smi)
            results.append({"index": index, "affinity": None,
                            "status": "invalid_input" if error else "ok", "error": error})
            if error is None:
                valid.append((index, mol))
        for start in range(0, len(valid), self.batch_size):
            batch = valid[start:start + self.batch_size]
            values = self.predictor.predict([mol for _, mol in batch])
            if len(values) != len(batch):
                raise RuntimeError("Predictor returned the wrong number of affinities")
            for (index, _), value in zip(batch, values):
                results[index]["affinity"] = finite_number(value, "predicted affinity")
        return {"version": PROTOCOL_VERSION, "target_id": self.target_id,
                "oracle_id": self.predictor.oracle_id, "results": results}


def make_server(oracle, port=8765):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):
            LOGGER.info(fmt, *args)

        def send_json(self, code, data):
            body = json.dumps(data, allow_nan=False).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                LOGGER.warning("Client disconnected before response")

        def do_GET(self):
            if self.path != "/health":
                return self.send_json(404, {"error": "not_found"})
            self.send_json(200, oracle.health())

        def do_POST(self):
            if self.path != "/score":
                return self.send_json(404, {"error": "not_found"})
            try:
                size = int(self.headers.get("Content-Length", "-1"))
                if size < 0 or self.headers.get("Transfer-Encoding"):
                    raise ValueError("Content-Length required; chunked requests unsupported")
                data = loads(self.rfile.read(size))
                if not isinstance(data, dict) or set(data) != {"version", "target_id", "smiles"}:
                    raise ValueError("Expected version, target_id and smiles")
                if type(data["version"]) is not int or data["version"] != PROTOCOL_VERSION:
                    raise ValueError("Unsupported protocol version")
                if not isinstance(data["target_id"], str) or not data["target_id"].strip():
                    raise ValueError("target_id must be a nonempty string")
                if not isinstance(data["smiles"], list) or not all(isinstance(s, str) for s in data["smiles"]):
                    raise ValueError("smiles must be an array of strings")
            except (ValueError, UnicodeError) as exc:
                return self.send_json(400, {"error": str(exc)})
            if data["target_id"] != oracle.target_id:
                return self.send_json(409, {"error": "target_id_mismatch"})
            try:
                response = oracle.score(data["smiles"])
            except Exception:
                LOGGER.exception("PLANET inference failed")
                return self.send_json(500, {"error": "inference_failed"})
            self.send_json(200, response)

    # HTTPServer intentionally processes one request at a time, including model inference.
    return HTTPServer(("127.0.0.1", port), Handler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--manifest", help="New output JSON file; refuses to overwrite")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        from .backend import PlanetPredictor

        cfg = server_config(args.config)
        predictor = PlanetPredictor(cfg)
        manifest = json.dumps(predictor.manifest, indent=2, sort_keys=True, allow_nan=False)
        if args.manifest:
            with Path(args.manifest).open("x", encoding="utf-8") as handle:
                handle.write(manifest + "\n")
        LOGGER.info("Resolved oracle manifest: %s", manifest)
        with make_server(Oracle(predictor, cfg["target_id"], cfg["batch_size"]), cfg["port"]) as server:
            LOGGER.info("Ready at http://127.0.0.1:%s", cfg["port"])
            try:
                server.serve_forever()
            except KeyboardInterrupt:
                LOGGER.info("Stopping oracle")
    except Exception:
        LOGGER.exception("Oracle startup failed")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
