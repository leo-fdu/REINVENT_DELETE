import threading

import pytest
from rdkit import Chem

from planet_oracle.chemistry import prepare_molecule
from planet_oracle.server import Oracle, make_server

class FakePredictor:
    oracle_id = "test-oracle"
    device = "cpu"
    pocket_info = {"residue_count": 1}

    def __init__(self):
        self.batches = []

    def prepare(self, smiles):
        return prepare_molecule(smiles, max_neighbors=10)

    def predict(self, molecules):
        self.batches.append([Chem.MolToSmiles(Chem.RemoveHs(mol)) for mol in molecules])
        return [float(4 + mol.GetNumHeavyAtoms()) for mol in molecules]


@pytest.fixture
def running_oracle():
    predictor = FakePredictor()
    oracle = Oracle(predictor, "test-target", batch_size=2)
    with make_server(oracle, port=0) as server:
        thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
        thread.start()
        cfg = {"url": f"http://127.0.0.1:{server.server_port}", "target_id": "test-target",
               "oracle_id": predictor.oracle_id, "timeout": 2}
        try:
            yield oracle, cfg
        finally:
            server.shutdown()
            thread.join(timeout=5)
