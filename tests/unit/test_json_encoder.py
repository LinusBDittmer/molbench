import json

import numpy
import pytest

from molbench.json_encoder import MolbenchJSONEncoder
from molbench.molecule import Datapoint, Molecule


def _dumps(obj):
    return json.dumps(obj, cls=MolbenchJSONEncoder)


def test_encodes_datapoint():
    dp = Datapoint(1.5, "eV")
    result = json.loads(_dumps(dp))
    assert result == {"value": 1.5, "unit": "eV"}


def test_encodes_datapoint_negative():
    dp = Datapoint(-76.0, "au")
    result = json.loads(_dumps(dp))
    assert result["value"] == pytest.approx(-76.0)
    assert result["unit"] == "au"


def test_encodes_numpy_array():
    arr = numpy.array([1.0, 2.0, 3.0])
    result = json.loads(_dumps(arr))
    assert result == [1.0, 2.0, 3.0]


def test_encodes_numpy_int():
    val = numpy.int64(42)
    result = json.loads(_dumps(val))
    assert result == 42


def test_encodes_plain_int():
    result = json.loads(_dumps(42))
    assert result == 42


def test_encodes_plain_str():
    result = json.loads(_dumps("hello"))
    assert result == "hello"


def test_encodes_molecule():
    mol = Molecule(
        "water",
        "bench",
        {"xyz": "O 0 0 0", "charge": 0},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": Datapoint(-76.0, "au")},
            }
        },
    )
    result = json.loads(_dumps(mol))
    assert "water" in result
    assert result["water"]["charge"] == 0


def test_datapoint_roundtrip():
    dp = Datapoint(3.14, "eV")
    serialized = _dumps(dp)
    restored = json.loads(serialized)
    assert Datapoint(restored["value"], restored["unit"]) == dp


def test_list_of_datapoints():
    data = [Datapoint(1.0, "au"), Datapoint(2.0, "au")]
    result = json.loads(_dumps(data))
    assert result[0] == {"value": 1.0, "unit": "au"}
    assert result[1] == {"value": 2.0, "unit": "au"}


def test_encoding_molecule_does_not_mutate_system_data():
    mol = Molecule(
        "water",
        "bench",
        {"xyz": "O 0 0 0", "charge": 0},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": Datapoint(-76.0, "au")},
            }
        },
    )
    before = dict(mol.system_data)
    _dumps(mol)
    assert mol.system_data == before
    assert "properties" not in mol.system_data
