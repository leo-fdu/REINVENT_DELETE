#!/usr/bin/env python3
"""Local browser interface for recording manual ligand task designs."""

from __future__ import annotations

from argparse import ArgumentParser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit
import json
import re

import design

HERE = Path(__file__).resolve().parent
STATIC = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/style.css": ("style.css", "text/css; charset=utf-8"),
}
RECORD_PATH = re.compile(r"/(20\d\d-\d\d-\d\dT\d\d-\d\d-\d\d-\d{6}[+-]\d{4})/(index\.html|designs\.json|figures/[a-z0-9]+_(?:libinvent|linkinvent)\.svg)")


class Handler(BaseHTTPRequestHandler):
    def send_bytes(self, status: int, payload: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def send_json(self, status: int, value: object) -> None:
        self.send_bytes(status, json.dumps(value, ensure_ascii=False).encode("utf-8"), "application/json; charset=utf-8")

    def do_GET(self) -> None:
        path = unquote(urlsplit(self.path).path)
        if path in STATIC:
            filename, content_type = STATIC[path]
            self.send_bytes(200, (HERE / filename).read_bytes(), content_type)
            return
        if path == "/api/targets":
            self.send_json(200, [
                {"target": target, "pdb": design.atlas.PDB[target], "sha256": design.atlas.RAW_SHA256[target]}
                for target in design.atlas.CODES
            ])
            return
        if path.startswith("/api/targets/"):
            try:
                self.send_json(200, design.target_info(path.removeprefix("/api/targets/")))
            except ValueError as error:
                self.send_json(404, {"error": str(error)})
            return
        match = RECORD_PATH.fullmatch(path)
        if match:
            folder, filename = match.groups()
            record_file = HERE / folder / filename
            if record_file.is_file():
                content_type = ("text/html" if filename.endswith(".html") else
                                "application/json" if filename.endswith(".json") else "image/svg+xml")
                self.send_bytes(200, record_file.read_bytes(), content_type + "; charset=utf-8")
                return
        self.send_json(404, {"error": "页面不存在"})

    def do_POST(self) -> None:
        path = urlsplit(self.path).path
        if path not in ("/api/preview", "/api/record"):
            self.send_json(404, {"error": "接口不存在"})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= 1_000_000:
                raise ValueError("请求体大小无效")
            body = json.loads(self.rfile.read(length))
            if not isinstance(body, dict):
                raise ValueError("请求体必须是 JSON 对象")
            if path == "/api/preview":
                result = design.preview(body.get("target"), body.get("mode"), body.get("cuts"), body.get("retained"))
                self.send_json(200, result)
            else:
                folder = design.save_snapshot(body.get("designs"))
                self.send_json(201, {"directory": folder.name, "report": f"/{folder.name}/index.html"})
        except (ValueError, AssertionError, KeyError, TypeError, json.JSONDecodeError) as error:
            self.send_json(400, {"error": str(error)})
        except OSError as error:
            self.send_json(500, {"error": f"写入记录失败：{error}"})


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8766, help="本机端口，默认 8766")
    args = parser.parse_args()
    for target in design.atlas.CODES:
        design.ligand(target)
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"手动任务设计界面：http://127.0.0.1:{args.port}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
