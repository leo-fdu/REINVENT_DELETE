import io
import json
import subprocess
import sys
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from planet_oracle.client import (
    INVALID_AFFINITY_SENTINEL,
    client_config,
    score,
    stdin_smiles,
    validate_response,
)
from planet_oracle.config import loads, server_config


def run_client(tmp_path, cfg, stdin):
    config = tmp_path / "client.json"
    config.write_text(json.dumps(cfg))
    return subprocess.run([sys.executable, "-S", "-m", "planet_oracle.client", "--config", str(config)],
                          input=stdin, text=True, capture_output=True, timeout=10)


def test_batches_order_duplicates_and_empty(running_oracle):
    oracle, cfg = running_oracle
    result = score(["CCO", "invalid", "", "C*", "C.C", "C", "CCO", "CC"], cfg)["payload"]
    assert result["planet_affinity"] == [7., None, None, None, None, 5., 7., 6.]
    assert result["planet_status"] == ["ok"] + ["invalid_input"] * 4 + ["ok"] * 3
    assert result["planet_affinity_for_scoring"] == [
        7., *([INVALID_AFFINITY_SENTINEL] * 4), 5., 7., 6.
    ]
    assert all(len(values) == 8 for values in result.values())
    assert oracle.predictor.batches == [["CCO", "C"], ["CCO", "CC"]]
    assert all(values == [] for values in score([], cfg)["payload"].values())
    assert score(["invalid"], cfg)["payload"]["planet_affinity"] == [None]
    assert len(oracle.predictor.batches) == 2
    with urlopen(cfg["url"] + "/health") as response:
        assert json.load(response) == oracle.health()


@pytest.mark.parametrize("text,expected", [("", []), ("\n", [""]), ("CCO\n\nC\n", ["CCO", "", "C"]),
                                         ("CCO\r\nC", ["CCO", "C"])])
def test_stdin(text, expected):
    assert stdin_smiles(io.StringIO(text)) == expected


def test_client_subprocess_standard_library_only(running_oracle, tmp_path):
    _, cfg = running_oracle
    result = run_client(tmp_path, cfg, "CCO\n\nCCO\n")
    assert result.returncode == 0, result.stderr
    assert result.stderr == ""
    assert json.loads(result.stdout)["payload"]["planet_affinity"] == [7., None, 7.]
    result = run_client(tmp_path, cfg, "")
    assert result.returncode == 0
    assert json.loads(result.stdout)["payload"]["planet_affinity_for_scoring"] == []


@pytest.mark.parametrize("change", [{"target_id": "wrong"}, {"oracle_id": "wrong"},
                                     {"url": "http://127.0.0.1:0"}])
def test_client_failures_nonzero(running_oracle, tmp_path, change):
    _, cfg = running_oracle
    result = run_client(tmp_path, {**cfg, **change}, "CCO")
    assert result.returncode != 0
    assert result.stdout == ""
    assert "PLANET client failed" in result.stderr


@pytest.mark.parametrize("body,status", [
    (b"{", 400), (b"[]", 400),
    (b'{"version":true,"target_id":"test-target","smiles":[]}', 400),
    (b'{"version":2,"target_id":"test-target","smiles":[]}', 400),
    (b'{"version":1,"target_id":"test-target","smiles":[null]}', 400),
    (b'{"version":1,"target_id":"test-target","smiles":[NaN]}', 400),
    (b'{"version":1,"target_id":"test-target","smiles":"CCO"}', 400),
    (b'{"version":1,"smiles":[]}', 400),
    (b'{"version":1,"target_id":"wrong","smiles":[]}', 409),
])
def test_http_validation(running_oracle, body, status):
    _, cfg = running_oracle
    with pytest.raises(HTTPError) as error:
        urlopen(Request(cfg["url"] + "/score", data=body))
    assert error.value.code == status


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf"), True, "7"])
def test_bad_prediction_fails_whole_request(running_oracle, monkeypatch, tmp_path, bad):
    oracle, cfg = running_oracle
    monkeypatch.setattr(oracle.predictor, "predict", lambda mols: [bad] * len(mols))
    result = run_client(tmp_path, cfg, "CCO\ninvalid\nCC")
    assert result.returncode != 0 and result.stdout == ""
    assert "500" in result.stderr


@pytest.mark.parametrize("failure", [RuntimeError("model failed"), MemoryError("out of memory")])
def test_model_failure(running_oracle, monkeypatch, failure):
    oracle, cfg = running_oracle

    def fail(molecules):
        raise failure

    monkeypatch.setattr(oracle.predictor, "predict", fail)
    with pytest.raises(HTTPError) as error:
        score(["CCO"], cfg)
    assert error.value.code == 500


def test_later_batch_failure_and_wrong_length(running_oracle, monkeypatch):
    oracle, cfg = running_oracle
    original = oracle.predictor.predict

    def fail_second(molecules):
        if oracle.predictor.batches:
            raise RuntimeError("batch graph construction failed")
        return original(molecules)

    monkeypatch.setattr(oracle.predictor, "predict", fail_second)
    with pytest.raises(HTTPError) as error:
        score(["CCO", "CC", "C"], cfg)
    assert error.value.code == 500
    monkeypatch.setattr(oracle.predictor, "predict", lambda molecules: [])
    with pytest.raises(HTTPError) as error:
        score(["CCO"], cfg)
    assert error.value.code == 500


def test_invalid_response_json_over_http(running_oracle, monkeypatch, tmp_path):
    # Patch only the server handler for this independent server, not global json.dumps.
    from planet_oracle.server import make_server
    import threading

    oracle, cfg = running_oracle
    with make_server(oracle, port=0) as server:
        def bad_response(handler):
            handler.send_response(200)
            handler.send_header("Content-Length", "1")
            handler.end_headers()
            handler.wfile.write(b"{")

        monkeypatch.setattr(server.RequestHandlerClass, "do_POST", bad_response)
        thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
        thread.start()
        try:
            result = run_client(tmp_path, {**cfg, "url": f"http://127.0.0.1:{server.server_port}"}, "CCO")
            assert result.returncode != 0 and result.stdout == ""
            assert "JSONDecodeError" in result.stderr
        finally:
            server.shutdown()
            thread.join(timeout=5)


def test_timeout(running_oracle, monkeypatch, tmp_path):
    oracle, cfg = running_oracle
    original = oracle.predictor.predict

    def slow(molecules):
        time.sleep(0.15)
        return original(molecules)

    monkeypatch.setattr(oracle.predictor, "predict", slow)
    result = run_client(tmp_path, {**cfg, "timeout": 0.02}, "CCO")
    assert result.returncode != 0 and result.stdout == ""
    assert "timed out" in result.stderr


@pytest.mark.parametrize("mutation", [
    lambda d: d.update(version=True), lambda d: d.update(version=2),
    lambda d: d.update(target_id="wrong"), lambda d: d.update(oracle_id="wrong"),
    lambda d: d.update(results=[]), lambda d: d["results"].append(d["results"][0]),
    lambda d: d["results"][0].pop("index"),
    lambda d: d["results"][0].update(index=True),
    lambda d: d["results"][1].update(index=0),
    lambda d: d["results"].reverse(),
    lambda d: d["results"][0].update(affinity=float("nan")),
    lambda d: d["results"][0].update(affinity=True),
    lambda d: d["results"][0].update(error="failed"),
    lambda d: d["results"][0].update(status="failed"),
    lambda d: d["results"][0].update(status="invalid_input", error="invalid_smiles"),
    lambda d: d["results"][0].update(status="invalid_input", affinity=None),
])
def test_response_validation(mutation):
    cfg = {"target_id": "test-target", "oracle_id": "test-oracle"}
    data = {"version": 1, **cfg, "results": [
        {"index": i, "affinity": 7., "status": "ok", "error": None} for i in range(2)]}
    assert len(validate_response(data, 2, cfg)) == 2
    mutation(data)
    with pytest.raises(ValueError):
        validate_response(data, 2, cfg)


def test_malformed_response_is_nonzero(running_oracle, monkeypatch, tmp_path):
    oracle, cfg = running_oracle
    monkeypatch.setattr(oracle, "score", lambda smiles: {"version": 1, "results": []})
    result = run_client(tmp_path, cfg, "CCO")
    assert result.returncode != 0 and result.stdout == ""


@pytest.mark.parametrize("text", ["{", '{"x":NaN}', '{"x":Infinity}'])
def test_bad_json(text):
    with pytest.raises(ValueError):
        loads(text)


def test_server_config_resolution(tmp_path):
    for name in ("PLANET.param", "protein.pdb", "ligand.sdf"):
        (tmp_path / name).touch()
    cfg = {"target_id": "t", "planet_root": ".", "protein_pdb": "protein.pdb", "ligand_sdf": "ligand.sdf"}
    path = tmp_path / "server.json"
    path.write_text(json.dumps(cfg))
    resolved = server_config(path)
    assert resolved["checkpoint"] == str(tmp_path / "PLANET.param")
    assert resolved["protein_pdb"] == str(tmp_path / "protein.pdb")
    assert (resolved["device"], resolved["seed"], resolved["batch_size"]) == ("cpu", 42, 32)
    for update in ({"center": [0, 0, 0]}, {"batch_size": 0}, {"port": True}, {"unknown": 1},
                   {"seed": -1}, {"protein_pdb": "missing"}):
        path.write_text(json.dumps({**cfg, **update}))
        with pytest.raises(ValueError):
            server_config(path)
    del cfg["ligand_sdf"]
    for center in (None, [], [0, 1], [0, 0, True], [0, 0, float("inf")]):
        path.write_text(json.dumps({**cfg, "center": center}))
        with pytest.raises(ValueError):
            server_config(path)
    path.write_text(json.dumps({**cfg, "center": [0, 1, 2]}))
    assert server_config(path)["center"] == [0, 1, 2]


def test_client_config_defaults_and_validation(tmp_path):
    path = tmp_path / "client.json"
    cfg = {"url": "http://127.0.0.1:8765/", "target_id": "t", "oracle_id": "abc"}
    path.write_text(json.dumps(cfg))
    assert client_config(path)["timeout"] == 120
    for update in ({"url": "file:///tmp/file"}, {"oracle_id": ""}, {"target_id": None},
                   {"timeout": 0}, {"low": 4}, {"unknown": 2}):
        path.write_text(json.dumps({**cfg, **update}))
        with pytest.raises(ValueError):
            client_config(path)
