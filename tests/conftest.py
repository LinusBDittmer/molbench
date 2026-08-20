import json

import pytest

from molbench.comparison import Comparison
from molbench.molecule import Datapoint, Molecule, MoleculeList

# ---------------------------------------------------------------------------
# Minimal molecule fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def hydrogen_molecule():
    return Molecule(
        name="H",
        data_id="ref",
        system_data={"xyz": "H 0.0 0.0 0.0", "charge": 0, "multiplicity": 2},
        state_data={
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": Datapoint(-0.5, "au")},
            }
        },
    )


@pytest.fixture
def water_molecule():
    return Molecule(
        name="water",
        data_id="ref",
        system_data={
            "xyz": "O 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0",
            "charge": 0,
            "multiplicity": 1,
        },
        state_data={
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": Datapoint(-76.0, "au")},
            }
        },
    )


@pytest.fixture
def simple_molecule_list(hydrogen_molecule, water_molecule):
    ml = MoleculeList()
    ml.extend([hydrogen_molecule, water_molecule])
    return ml


@pytest.fixture
def simple_comparison(simple_molecule_list):
    c = Comparison()
    c.add(simple_molecule_list)
    return c


# ---------------------------------------------------------------------------
# Benchmark JSON fixtures on disk
# ---------------------------------------------------------------------------

MINIMAL_BENCHMARK = {
    "molA": {
        "charge": 0,
        "multiplicity": 1,
        "n_atoms": 1,
        "xyz": ["H 0.0 0.0 0.0"],
        "properties": {
            "s1": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": {"value": -0.5, "unit": "au"}},
                "transition_id": "s0->s1",
            }
        },
    },
    "molB": {
        "charge": 0,
        "multiplicity": 1,
        "n_atoms": 3,
        "xyz": ["O 0 0 0", "H 0 0 1", "H 0 1 0"],
        "properties": {
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": {"value": -76.0, "unit": "au"}},
            }
        },
    },
}

MINIMAL_BENCHMARK_LIST = {
    "relative": {
        "charge_list": [0, 1],
        "multiplicity_list": [1, 2],
        "n_atoms_list": [1, 1],
        "xyz_list": [["H 0 0 0"], ["H 0 0 1"]],
        "properties": {
            "p1": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": {"value": 0.5, "unit": "au"}},
                "stochiometry": [-1.0, 1.0],
            }
        },
    }
}


@pytest.fixture
def minimal_benchmark_file(tmp_path):
    f = tmp_path / "bench.json"
    f.write_text(json.dumps(MINIMAL_BENCHMARK))
    return str(f)


@pytest.fixture
def minimal_benchmark_list_file(tmp_path):
    f = tmp_path / "bench_list.json"
    f.write_text(json.dumps(MINIMAL_BENCHMARK_LIST))
    return str(f)


@pytest.fixture
def simple_template_file(tmp_path):
    content = "charge=[[charge]]\nbasis=[[basis]]\nxyz=\n[[xyz]]"
    f = tmp_path / "tmpl.txt"
    f.write_text(content)
    return str(f)


# ---------------------------------------------------------------------------
# Known-delta Comparison fixture for statistics tests
# ---------------------------------------------------------------------------


@pytest.fixture
def known_comparison():
    """Comparison with ref=-76.0 au (TBE) and interest=-75.9 au (HF).
    Signed error = interest - reference = +0.1 au.
    """
    ref_mol = Molecule(
        "water",
        "ref",
        {},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": Datapoint(-76.0, "au")},
            }
        },
    )
    int_mol = Molecule(
        "water",
        "computed",
        {},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": Datapoint(-75.9, "au")},
            }
        },
    )
    c = Comparison()
    c.add(MoleculeList([ref_mol, int_mol]))
    return c


@pytest.fixture
def two_molecule_comparison():
    """Two molecules, errors +0.1 and +0.3 au.
    MSE=0.2, MAE=0.2, RMSD=sqrt(0.05)≈0.2236, min=0.1, max=0.3
    """
    mols = MoleculeList(
        [
            Molecule(
                "water",
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(-76.0, "au")},
                    }
                },
            ),
            Molecule(
                "benzene",
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(-10.0, "au")},
                    }
                },
            ),
            Molecule(
                "water",
                "computed",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-75.9, "au")},
                    }
                },
            ),
            Molecule(
                "benzene",
                "computed",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-9.7, "au")},
                    }
                },
            ),
        ]
    )
    c = Comparison()
    c.add(mols)
    return c
